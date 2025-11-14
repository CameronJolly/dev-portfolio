struct PressureGlobals {
    gasConst: f32,
    restDensity: f32,
    viscosityCoeff: f32,
    h: f32,
    h2: f32,
    h5: f32,
    h8: f32,
    numParticles: u32,
};

struct Bounds{
    minX: f32,
    minY: f32,
    maxX: f32,
    maxY: f32,
    width: f32,
    height: f32,
    padding: vec2<f32>
}

@group(0) @binding(0) var<uniform> bounds: Bounds;
@group(0) @binding(1) var<uniform> globals: PressureGlobals;
@group(0) @binding(2) var<storage, read_write> densities: array<f32>;
@group(0) @binding(3) var<storage, read_write> positions: array<f32>;
@group(0) @binding(4) var<storage, read_write> accelerations: array<f32>;
@group(0) @binding(5) var<storage, read_write> viscocities: array<f32>;
@group(0) @binding(6) var<storage, read_write> velocities: array<f32>;

const PARTICLE_RADIUS: f32 = 0.03;
const RESTITUTION: f32 = 0.5;
const WALL_FRICTION: f32 = 0.1;
const DENSITY_EPS: f32 = 1e-6;

fn poly6Kernel(r2: f32, h2: f32) -> f32 {
    let constant = 4.0 / (3.14159265359 * globals.h8);
    if (r2 >= 0.0 && r2 <= h2) {
        let term = h2 - r2;
        return constant * term * term * term;
    }
    return 0.0;
}

fn spikyKernel(r2: f32, h: f32) -> f32 {
    let dist = sqrt(r2);

    if (dist == 0.0 || dist >= h){
        return 0.0;
    }

    let constant = -10.0 / (3.14159265359 * globals.h5);
    let term = (h - dist) * (h - dist);
    return constant * term / dist;
}

fn getCellNum(posX: f32, posY: f32) -> f32 {
    let cellX = floor(((posX*2.0) - bounds.maxX) / globals.h);
    let cellY = floor((((posY*2.0)+1.0) - bounds.maxY) / globals.h);
    return cellY*bounds.width + cellX;
}

fn calculateDensity(global_id: vec3<u32>){
    let epislon = 1e-4;
        
    const ix = positions[global_id.x * 2];
    const iy = positions[global_id.x * 2 + 1];

    let rho = 0.0;

    for(let index: u32 = 0; index < globals.numParticles; index += 1){
        if(index == global_id.x){continue;}

        const jx = positions[index * 2];
        const jy = positions[index * 2 + 1];

        const dx = ix-jx;
        const dy = iy-jy;

        const r2 = dx*dx + dy*dy;

        if (index != global_id.x && r2 < epislon * epislon) {
            const angle = 0.1 * 3.14159265359 * 2.0;
            
            positions[global_id.x * 2] += cos(angle) * epislon * 0.5;
            positions[global_id.x * 2 + 1] += sin(angle) * epislon * 0.5;
            positions[index * 2] -= cos(angle) * epislon * 0.5;
            positions[index * 2 + 1] -= sin(angle) * epislon * 0.5;
        }

        if(r2 < globals.h2){
            rho += poly6Kernel(r2, globals.h2);
        }

    }
    workgroupBarrier();
    densities[global_id.x] += rho;
}

fn calculatePressure(global_id: vec3<u32>){
    if (global_id.x >= arrayLength(&resultArray)) {
        return;
    }
    let density = densities[global_id.x];
    let pressure = max(0.0, globals.gasConst * (density - globals.restDensity));
    resultArray[global_id.x] = pressure;
}

fn calculatePressureForceAndViscocity(global_id: vec3<u32>){
    if (global_id.x >= globals.numParticles) {
        return;
    }

    let idx = global_id.x;
    let ix = positions[idx * 2];
    let iy = positions[idx * 2 + 1];
    let ivx = velocities[idx * 2];
    let ivy = velocities[idx * 2 + 1];
    let density_i = max(densities[idx], DENSITY_EPS);
    let pressure_i = resultArray[idx];

    var fx = 0.0;
    var fy = 0.0;
    var viscX = 0.0;
    var viscY = 0.0;

    for (var j: u32 = 0u; j < globals.numParticles; j = j + 1u) {
        if (j == idx) {
            continue;
        }

        let jx = positions[j * 2];
        let jy = positions[j * 2 + 1];
        let dx = ix - jx;
        let dy = iy - jy;
        let r2 = dx * dx + dy * dy;

        if (r2 > 0.0 && r2 <= globals.h2) {
            let gradMagnitude = spikyKernel(r2, globals.h);
            let grad = vec2<f32>(dx * gradMagnitude, dy * gradMagnitude);

            let density_j = max(densities[j], DENSITY_EPS);
            let pressure_j = resultArray[j];
            let coeff = (pressure_i / (density_i * density_i)) + (pressure_j / (density_j * density_j));

            fx -= coeff * grad.x;
            fy -= coeff * grad.y;

            let jvx = velocities[j * 2];
            let jvy = velocities[j * 2 + 1];
            let w = poly6Kernel(r2, globals.h2);
            let invDensity_j = 1.0 / density_j;
            let velDiff = vec2<f32>(jvx - ivx, jvy - ivy);

            viscX += velDiff.x * invDensity_j * w;
            viscY += velDiff.y * invDensity_j * w;
        }
    }

    accelerations[idx * 2] += fx + (viscosityCoeff * viscX);
    accelerations[idx * 2 + 1] += fy + (viscosityCoeff * viscY);
}

fn handleCollision(global_id: vec3<u32>){
    if (global_id.x >= globals.numParticles) {
        return;
    }

    let idx = global_id.x;
    var px = positions[idx * 2];
    var py = positions[idx * 2 + 1];
    var vx = velocities[idx * 2];
    var vy = velocities[idx * 2 + 1];

    let minX = bounds.minX + PARTICLE_RADIUS;
    let maxX = bounds.maxX - PARTICLE_RADIUS;
    let minY = bounds.minY + PARTICLE_RADIUS;
    let maxY = bounds.maxY - PARTICLE_RADIUS;

    if (px < minX) {
        px = minX;
        vx = -vx * RESTITUTION;
        vy *= WALL_FRICTION;
    } else if (px > maxX) {
        px = maxX;
        vx = -vx * RESTITUTION;
        vy *= WALL_FRICTION;
    }

    if (py < minY) {
        py = minY;
        vy = -vy * RESTITUTION;
        vx *= WALL_FRICTION;
    } else if (py > maxY) {
        py = maxY;
        vy = -vy * RESTITUTION;
        vx *= WALL_FRICTION;
    }

    positions[idx * 2] = px;
    positions[idx * 2 + 1] = py;
    velocities[idx * 2] = vx;
    velocities[idx * 2 + 1] = vy;
}


@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    if (global_id.x >= globals.numParticles) {
        return;
    }
    calculateDensity(global_id);
    workgroupBarrier();
    calculatePressure(global_id);
    workgroupBarrier();
    calculatePressureForceAndViscocity(global_id);
    handleCollision(global_id);
}