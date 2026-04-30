struct PressureGlobals {
    gasConst: f32,
    restDensity: f32,
    viscosityCoeff: f32,
    h: f32,
    h2: f32,
    h5: f32,
    h8: f32,
    numParticles: u32,
    gravity: f32,
    dt: f32,
    padding: vec2<f32>
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

@compute @workgroup_size(64)
fn densityPass(@builtin(global_invocation_id) global_id: vec3<u32>) {
    if (global_id.x >= globals.numParticles) {
        return;
    }

    let idx = global_id.x;
    let ix = positions[idx * 2u];
    let iy = positions[idx * 2u + 1u];
    let epsilon = 1e-4;

    var rho = 0.0;

    for (var j: u32 = 0u; j < globals.numParticles; j += 1u) {
        if (j == idx) { continue; }

        let jBase = j * 2u;
        let jx = positions[jBase];
        let jy = positions[jBase + 1u];

        let dx = ix - jx;
        let dy = iy - jy;
        let r2 = dx * dx + dy * dy;

        if (r2 < globals.h2) {
            rho += poly6Kernel(r2, globals.h2);
        }
    }
    densities[idx] = rho;
}

@compute @workgroup_size(64)
fn forcePass(@builtin(global_invocation_id) global_id: vec3<u32>) {
    if (global_id.x >= globals.numParticles) {
        return;
    }

    let idx = global_id.x;
    let ix = positions[idx * 2u];
    let iy = positions[idx * 2u + 1u];
    let ivx = velocities[idx * 2u];
    let ivy = velocities[idx * 2u + 1u];
    let density_i = max(densities[idx], DENSITY_EPS);
    let pressure_i = max(0.0, globals.gasConst * (density_i - globals.restDensity));

    var fx = 0.0;
    var fy = globals.gravity;
    var viscX = 0.0;
    var viscY = 0.0;

    for (var j: u32 = 0u; j < globals.numParticles; j += 1u) {
        if (j == idx) {
            continue;
        }

        let jBase = j * 2u;
        let jx = positions[jBase];
        let jy = positions[jBase + 1u];
        let dx = ix - jx;
        let dy = iy - jy;
        let r2 = dx * dx + dy * dy;

        if (r2 > 0.0 && r2 <= globals.h2) {
            let gradMagnitude = spikyKernel(r2, globals.h);
            let grad = vec2<f32>(dx * gradMagnitude, dy * gradMagnitude);

            let density_j = max(densities[j], DENSITY_EPS);
            let pressure_j = max(0.0, globals.gasConst * (density_j - globals.restDensity));
            let coeff = (pressure_i / (density_i * density_i)) + (pressure_j / (density_j * density_j));

            fx -= coeff * grad.x;
            fy -= coeff * grad.y;

            let jvx = velocities[j * 2u];
            let jvy = velocities[j * 2u + 1u];
            let w = poly6Kernel(r2, globals.h2);
            let invDensity_j = 1.0 / density_j;

            viscX += (jvx - ivx) * invDensity_j * w;
            viscY += (jvy - ivy) * invDensity_j * w;
        }
    }

    accelerations[idx * 2u] = fx + (globals.viscosityCoeff * viscX);
    accelerations[idx * 2u + 1u] = fy + (globals.viscosityCoeff * viscY);
}

@compute @workgroup_size(64)
fn integratePass(@builtin(global_invocation_id) global_id: vec3<u32>) {
    if (global_id.x >= globals.numParticles) {
        return;
    }

    let idx = global_id.x;
    let dt = globals.dt;

    velocities[idx * 2u] += accelerations[idx * 2u] * dt;
    velocities[idx * 2u + 1u] += accelerations[idx * 2u + 1u] * dt;

    positions[idx * 2u] += velocities[idx * 2u] * dt;
    positions[idx * 2u + 1u] += velocities[idx * 2u + 1u] * dt;

    // Collision
    var px = positions[idx * 2u];
    var py = positions[idx * 2u + 1u];
    var vx = velocities[idx * 2u];
    var vy = velocities[idx * 2u + 1u];

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

    positions[idx * 2u] = px;
    positions[idx * 2u + 1u] = py;
    velocities[idx * 2u] = vx;
    velocities[idx * 2u + 1u] = vy;
}