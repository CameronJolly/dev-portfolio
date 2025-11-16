import * as THREE from 'three';

class particle extends THREE.Mesh {
    index;
    mass = 1;
    constructor(geometry, material) {
        super(geometry,material);
    }
}

class Engine {
    #particles = [];
    numParticles;

    forceDirection = true;

    positions;
    velocities;
    accelerations;
    pressures;
    viscocities;
    densities;

    cells = [];
    gridWidth;
    gridHeight;

    mouseX = 0;
    mouseY = 0;

    // Simulation constants
    H = .3;            // also grid cell size
    GAS_CONST = 80;
    GRAVITY_Y = -9;
    DENSITY_EPS = 1e-6;      // to avoid divide-by-zero
    PARTICLE_RADIUS = .03;

    RESTITUTION = 0.5;
    WALL_FRICTION = 0.1;
    REST_DENSITY = 4;
    VISCOSITY_COEFF = 80;

    // Fixed-step physics update for stability
    FIXED_DT = 1 / 100;   // seconds per physics step
    MAX_ACCUM = 0.05;     // clamp big frame gaps (s)
    MAX_STEPS = 8;        // prevent spiral of death

    lastTime = 0;
    accumulator = 0;

    //webgpu globals
    queue;
    device;
    densitiesBuffer;
    positionsBuffer;
    accelerationsBuffer;
    viscocitiesBuffer;
    velocitiesBuffer;
    sphCalcUniformBuffer;
    boundsUniformBuffer;
    sphCalcBindGroup;
    computePipeline;
    readbackBuffer;
    positionsBufferSize;
    positionsReadbackBuffer;



    async init() {
        await this.initWebGPU();
    }

    async initWebGPU() {
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported in this browser.");
        }

        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) {
            throw new Error("Couldn't request WebGPU adapter.");
        }

        this.device = await adapter.requestDevice();
        this.queue = this.device.queue;

        // getting wgsl text and parsing it for use in tsl pipeline
        const shaderUrl = new URL('./compute-shaders/sphPipeline.wgsl', import.meta.url);
        const shaderResponse = await fetch(shaderUrl);
        const shaderCode = await shaderResponse.text();

        // creating webgpu shader module which loads the raw wgsl an duses it
        const shaderModule = this.device.createShaderModule({code: shaderCode});

        //creating buffers to send to the gpu, they have to be strictly typed
        this.densitiesBuffer = this.device.createBuffer({
            size: Float32Array.BYTES_PER_ELEMENT * this.numParticles,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC
        });

        const positionsBytes = Float32Array.BYTES_PER_ELEMENT * this.numParticles * 2;
        this.positionsBufferSize = positionsBytes;
        this.positionsBuffer = this.device.createBuffer({
            size: positionsBytes,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC
        })

        this.accelerationsBuffer = this.device.createBuffer({
            size: Float32Array.BYTES_PER_ELEMENT * this.numParticles * 2,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST
        })

        this.viscocitiesBuffer = this.device.createBuffer({
            size: Float32Array.BYTES_PER_ELEMENT * this.numParticles,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST
        })

        
        this.velocitiesBuffer = this.device.createBuffer({
            size: Float32Array.BYTES_PER_ELEMENT * this.numParticles * 2,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST
        })

        //buffer for globals needed for the sph pipeline
        this.sphCalcUniformBuffer = this.device.createBuffer({
            size: 32,
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });

        this.boundsUniformBuffer = this.device.createBuffer({
            size: 32,
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });

        // creating layout for the pressure bind group (what is getting sent into the compute shader kernel)
        const pressureBindGroupLayout = this.device.createBindGroupLayout({
            entries: [
                {binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'uniform'}},
                {binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'uniform'}},
                {binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'storage'}},
                {binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'storage'}},
                {binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'storage'}},
                {binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'storage'}},
                {binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'storage'}},
            ]
        })

        this.sphCalcBindGroup = this.device.createBindGroup({
            layout: pressureBindGroupLayout,
            entries: [
                {binding: 0, resource: {buffer: this.boundsUniformBuffer}},
                {binding: 1, resource: {buffer: this.sphCalcUniformBuffer}},
                {binding: 2, resource: {buffer: this.densitiesBuffer}},
                {binding: 3, resource: {buffer: this.positionsBuffer}},
                {binding: 4, resource: {buffer: this.accelerationsBuffer}},
                {binding: 5, resource: {buffer: this.viscocitiesBuffer}},
                {binding: 6, resource: {buffer: this.velocitiesBuffer}},
            ]
        });

        //setup pressure compute pipeline
        this.computePipeline = this.device.createComputePipeline({
            layout: this.device.createPipelineLayout({
                bindGroupLayouts: [pressureBindGroupLayout]
            }),
            compute: {
                module: shaderModule,
                entryPoint: 'main'
            }
        });

        this.positionsReadbackBuffer = this.device.createBuffer({
            size: positionsBytes,
            usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
        });

        this.densitiesReadbackBuffer = this.device.createBuffer({
            size: Float32Array.BYTES_PER_ELEMENT * this.numParticles,
            usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
        });
    }

    poly6Kernel(r2, h) {
        const h2 = h * h;
        const h8 = h2 * h2 * h2 * h2;
        const constant = 4 / (Math.PI * h8);

        if (r2 >= 0 && r2 <= h2) {
            const term = h2 - r2;
            return constant * term * term * term;
        } else {
            return 0;
        }
    }

    spikyKernel(dx, dy, distSq, h) {
        const dist = Math.sqrt(distSq);

        if (dist === 0 || dist >= h) {
            return { gradWx: 0, gradWy: 0 };
        }

        const h5 = h * h * h * h * h;
        const constant = -10 / (Math.PI * h5);
        const term = (h - dist) * (h - dist);
        const gradMagnitude = constant * term / dist;

        return {
            gradWx: gradMagnitude * dx,
            gradWy: gradMagnitude * dy
        };
    }

    getNeighborCells(cellX, cellY) {
        const neighbors = [];
        for (let dy = -1; dy <= 1; dy++) {
            for (let dx = -1; dx <= 1; dx++) {
                const nx = cellX + dx;
                const ny = cellY + dy;
                if (nx >= 0 && nx < this.gridWidth && ny >= 0 && ny < this.gridHeight) {
                    const cellIndex = ny * this.gridWidth + nx;
                    neighbors.push(this.cells[cellIndex]);
                }
            }
        }
        return neighbors;
    }

    getNeighborCellsInRadius(cellX, cellY, radius) {
        const cellR = Math.ceil(radius / this.H);
        const minCX = Math.max(0, cellX - cellR);
        const maxCX = Math.min(this.gridWidth - 1, cellX + cellR);
        const minCY = Math.max(0, cellY - cellR);
        const maxCY = Math.min(this.gridHeight - 1, cellY + cellR);

        const cells = [];
        for (let y = minCY; y <= maxCY; y++) {
            for (let x = minCX; x <= maxCX; x++) {
                cells.push(this.cells[y * this.gridWidth + x]);
            }
        }
        return cells;
    }

    getCellCoords(particleIndex, bounds) {
        const px = this.positions[particleIndex * 2];
        const py = this.positions[particleIndex * 2 + 1];
        const cellX = Math.floor((px - bounds.minX) / this.H);
        const cellY = Math.floor((py - bounds.minY) / this.H);
        return { cellX, cellY };
    }

    createParticles(numParticles, minX, maxX, minY, maxY, scene) {
        // allocate gpu ready arrays for eventual sending into kernel
        this.#particles = [];

        this.positions = new Float32Array(numParticles*2);
        this.velocities = new Float32Array(numParticles*2);
        this.accelerations = new Float32Array(numParticles*2);

        this.pressures = new Float32Array(numParticles);
        this.viscocities = new Float32Array(numParticles);
        this.densities = new Float32Array(numParticles);

        // create grid
        this.gridWidth = Math.ceil((maxX - minX) / this.H);
        this.gridHeight = Math.ceil((maxY - minY) / this.H);
        const totalCells = this.gridWidth * this.gridHeight;

        for (let cellIndex = 0; cellIndex < totalCells ; cellIndex++) {
            this.cells.push([]);
        }

        //even distribution - grid layout
        const gridCols = Math.ceil(Math.sqrt(numParticles));
        const gridRows = Math.ceil(numParticles / gridCols);
        const spacingX = (maxX - minX) / (gridCols + 1);
        const spacingY = (maxY - minY) / (gridRows + 1);

        for (let index = 0; index < numParticles; index++) {
            const geometry = new THREE.CircleGeometry(0.07, 10);
            const material = new THREE.MeshBasicMaterial({ color: 0x4fd9f7 });
            let newParticle = new particle(geometry, material);

            const col = index % gridCols;
            const row = Math.floor(index / gridCols);

            newParticle.index = index;
            newParticle.position.x = minX + (col + 1) * spacingX;
            newParticle.position.y = minY + (row + 1) * spacingY;

            this.positions[index*2] = newParticle.position.x;
            this.positions[(index*2)+1] = newParticle.position.y;

            const cellY = Math.floor((newParticle.position.y-minY) / this.H);
            const cellX = Math.floor((newParticle.position.x-minX)/this.H);
            this.cells[(cellY*this.gridWidth)+cellX].push(newParticle.index);

            scene.add(newParticle);
            this.#particles.push(newParticle);
        }

        this.numParticles = this.#particles.length;
        return this.#particles;
    }

    applyPhysics(dt) {
        // Semi-implicit Euler
        this.#particles.forEach(p => {
            const idx = p.index;
            
            this.velocities[idx * 2] += this.accelerations[idx * 2] * dt;
            this.velocities[idx * 2 + 1] += this.accelerations[idx * 2 + 1] * dt;
            
            this.positions[idx * 2] += this.velocities[idx * 2] * dt;
            this.positions[idx * 2 + 1] += this.velocities[idx * 2 + 1] * dt;

            p.position.x = this.positions[idx * 2];
            p.position.y = this.positions[idx * 2 + 1];
        });
    }

    calculateDensity(bounds) {
        const h2 = this.H * this.H;
        const EPSILON = 1e-4; // Small separation distance
        
        this.#particles.forEach(p => {
            const targetElem = p.index;
            const { cellX, cellY } = this.getCellCoords(targetElem, bounds);
            const neighborCells = this.getNeighborCells(cellX, cellY);
            
            let rho = 0;
            neighborCells.forEach(cell => {
                cell.forEach(outerElem => {
                    const ox = this.positions[outerElem * 2];
                    const oy = this.positions[outerElem * 2 + 1];
                    const tx = this.positions[targetElem * 2];
                    const ty = this.positions[targetElem * 2 + 1];
                    const dx = ox - tx;
                    const dy = oy - ty;
                    const r2 = dx*dx + dy*dy;
                    
                    // Separate overlapping particles
                    if (outerElem !== targetElem && r2 < EPSILON * EPSILON) {
                        // Generate small random offset to prevent particles from staying stuck
                        const angle = Math.random() * Math.PI * 2;
                        const offsetDist = EPSILON;
                        
                        this.positions[targetElem * 2] += Math.cos(angle) * offsetDist * 0.5;
                        this.positions[targetElem * 2 + 1] += Math.sin(angle) * offsetDist * 0.5;
                        this.positions[outerElem * 2] -= Math.cos(angle) * offsetDist * 0.5;
                        this.positions[outerElem * 2 + 1] -= Math.sin(angle) * offsetDist * 0.5;
                    }
                    
                    if(r2 <= h2){
                        rho += this.poly6Kernel(r2, this.H);
                    }
                });
            });
            this.densities[targetElem] = rho;
        });
    }

    calculatePressure() {
        // this.queue.writeBuffer(this.densitiesBuffer, 0, this.densities);
        
        // const pressureCalcUniformData = new Float32Array([
        //     this.GAS_CONST,
        //     this.REST_DENSITY,
        //     0, 0
        // ]);
        // this.queue.writeBuffer(this.sphCalcUniformBuffer, 0, pressureCalcUniformData);

        // const pressureCommandEncoder = this.device.createCommandEncoder();
        // const passEncoder = pressureCommandEncoder.beginComputePass();

        // passEncoder.setPipeline(this.computePipeline);
        // passEncoder.setBindGroup(0, this.sphCalcBindGroup);

        // const workgroupCount = Math.ceil(this.numParticles / 64);
        // passEncoder.dispatchWorkgroups(workgroupCount);
        // passEncoder.end();

        // this.queue.submit([pressureCommandEncoder.finish()]);

        this.#particles.forEach(particle => {
            const pressure = this.GAS_CONST * (this.densities[particle.index]-this.REST_DENSITY);
            this.pressures[particle.index] = Math.max(0,pressure);
        });
        
    }

    calculatePressureForce(bounds) {
        const h2 = this.H * this.H;
        
        this.#particles.forEach(p => {
            const targetElem = p.index;
            const { cellX, cellY } = this.getCellCoords(targetElem, bounds);
            const neighborCells = this.getNeighborCells(cellX, cellY);
            
            let fx = 0, fy = 0;
            
            neighborCells.forEach(cell => {
                cell.forEach(outerElem => {
                    if (outerElem === targetElem) {
                        return;
                    }
                    
                    const ix = this.positions[targetElem * 2];
                    const iy = this.positions[targetElem * 2 + 1];
                    const jx = this.positions[outerElem * 2];
                    const jy = this.positions[outerElem * 2 + 1];
                    
                    const dx = ix - jx;
                    const dy = iy - jy;
                    const r2 = dx * dx + dy * dy;
                    
                    if (r2 > 0 && r2 <= h2) {
                        const { gradWx, gradWy } = this.spikyKernel(dx, dy, r2, this.H);
                        
                        const density_i = Math.max(this.densities[targetElem], this.DENSITY_EPS);
                        const density_j = Math.max(this.densities[outerElem], this.DENSITY_EPS);
                        const denom_i = density_i * density_i;
                        const denom_j = density_j * density_j;
                        
                        const pressure_i = this.pressures[targetElem];
                        const pressure_j = this.pressures[outerElem];
                        
                        const coeff = (pressure_i / denom_i + pressure_j / denom_j);
                        
                        fx -= coeff * gradWx;
                        fy -= coeff * gradWy;
                    }
                });
            });
            
            this.accelerations[targetElem * 2] += fx;
            this.accelerations[targetElem * 2 + 1] += fy;
        });
    }

    applyViscosity(dt, bounds) {
        const h2 = this.H * this.H;
        
        this.#particles.forEach(p => {
            const targetElem = p.index;
            const { cellX, cellY } = this.getCellCoords(targetElem, bounds);
            const neighborCells = this.getNeighborCells(cellX, cellY);
            
            let corrX = 0, corrY = 0;
            
            neighborCells.forEach(cell => {
                cell.forEach(outerElem => {
                    if (outerElem === targetElem) {
                        return;
                    }
                    
                    const ix = this.positions[targetElem * 2];
                    const iy = this.positions[targetElem * 2 + 1];
                    const jx = this.positions[outerElem * 2];
                    const jy = this.positions[outerElem * 2 + 1];
                    
                    const ivx = this.velocities[targetElem * 2];
                    const ivy = this.velocities[targetElem * 2 + 1];
                    const jvx = this.velocities[outerElem * 2];
                    const jvy = this.velocities[outerElem * 2 + 1];
                    
                    const dx = jx - ix;
                    const dy = jy - iy;
                    const r2 = dx * dx + dy * dy;
                    
                    if (r2 > 0 && r2 <= h2) {
                        const w = this.poly6Kernel(r2, this.H);
                        const densityJ = Math.max(this.densities[outerElem], this.DENSITY_EPS);
                        corrX += (jvx - ivx) * (1 / densityJ) * w;
                        corrY += (jvy - ivy) * (1 / densityJ) * w;
                    }
                });
            });
            
            this.velocities[targetElem * 2] += this.VISCOSITY_COEFF * corrX * dt;
            this.velocities[targetElem * 2 + 1] += this.VISCOSITY_COEFF * corrY * dt;
        });
    }

    setColorBasedOnDensity() {
        let minD = Infinity;
        let maxD = -Infinity;
        
        for (let i = 0; i < this.densities.length; i++) {
            if (this.densities[i] < minD) minD = this.densities[i];
            if (this.densities[i] > maxD) maxD = this.densities[i];
        }
        
        const range = Math.max(maxD - minD, 1e-8);

        this.#particles.forEach(p => {
            const t = (this.densities[p.index] - minD) / range;
            const red = t * 1.7;
            const green = 0;
            const blue = 1 - t;
            p.material.color.setRGB(red, green, blue);
        });
    }

    handleCollisions(bounds) {
        const minX = bounds.minX + this.PARTICLE_RADIUS;
        const maxX = bounds.maxX - this.PARTICLE_RADIUS;
        const minY = bounds.minY + this.PARTICLE_RADIUS;
        const maxY = bounds.maxY - this.PARTICLE_RADIUS;

        this.#particles.forEach(p => {
            const idx = p.index;
            let px = this.positions[idx * 2];
            let py = this.positions[idx * 2 + 1];
            let vx = this.velocities[idx * 2];
            let vy = this.velocities[idx * 2 + 1];

            // Left/right walls
            if (px < minX) {
                px = minX;
                vx = -vx * this.RESTITUTION;
                vy *= this.WALL_FRICTION;
            } else if (px > maxX) {
                px = maxX;
                vx = -vx * this.RESTITUTION;
                vy *= this.WALL_FRICTION;
            }

            // Bottom/top walls
            if (py < minY) {
                py = minY;
                vy = -vy * this.RESTITUTION;
                vx *= this.WALL_FRICTION;
            } else if (py > maxY) {
                py = maxY;
                vy = -vy * this.RESTITUTION;
                vx *= this.WALL_FRICTION;
            }

            // Write back to arrays and particle object
            this.positions[idx * 2] = px;
            this.positions[idx * 2 + 1] = py;
            this.velocities[idx * 2] = vx;
            this.velocities[idx * 2 + 1] = vy;
            
            p.position.x = px;
            p.position.y = py;
        });
    }

    rebuildGrid(bounds) {
        const minX = bounds.minX;
        const maxX = bounds.maxX;
        const minY = bounds.minY;
        const maxY = bounds.maxY;

        // Recalculate grid dimensions
        this.gridWidth = Math.ceil((maxX - minX) / this.H);
        this.gridHeight = Math.ceil((maxY - minY) / this.H);
        const totalCells = this.gridWidth * this.gridHeight;

        // Rebuild cell array
        this.cells = [];
        for (let cellIndex = 0; cellIndex < totalCells; cellIndex++) {
            this.cells.push([]);
        }
    }

    async stepPhysics(dt, bounds) {
        // Reset accelerations
        for (let i = 0; i < this.accelerations.length/2; i++) {
            this.accelerations[i*2] = 0;
            this.accelerations[(i*2)+1] = this.GRAVITY_Y;
        }

        this.rebuildGrid(bounds);
        this.cells.forEach(cell => cell.length = 0);
        this.#particles.forEach(p => {
            const { cellX, cellY } = this.getCellCoords(p.index, bounds);
            if (cellX >= 0 && cellX < this.gridWidth && cellY >= 0 && cellY < this.gridHeight) {
                this.cells[cellY * this.gridWidth + cellX].push(p.index);
            }
        });

        // SPH pipeline
        this.calculateDensity(bounds);
        this.calculatePressure();
        this.calculatePressureForce(bounds);
        this.applyViscosity(dt, bounds);
        this.handleCollisions(bounds);
        if(this.forceDirection){
            this.applyOutwardForce(bounds);
        }else{
            this.applyInwardForce(bounds);
        }

        this.applyPhysics(dt);
    }

    async applyWgslPhysics(dt, bounds){
        const h = this.H;
        const h2 = this.H * this.H;
        const h5 = h2 * h2 * this.H;
        const h8 = h5 * h2 * this.H;

        this.queue.writeBuffer(this.positionsBuffer, 0, this.positions);

        const globalsBuffer = new ArrayBuffer(32);
        const globalsView = new DataView(globalsBuffer);
        globalsView.setFloat32(0, this.GAS_CONST, true);
        globalsView.setFloat32(4, this.REST_DENSITY, true);
        globalsView.setFloat32(8, this.VISCOSITY_COEFF, true);
        globalsView.setFloat32(12, h, true);
        globalsView.setFloat32(16, h2, true);
        globalsView.setFloat32(20, h5, true);
        globalsView.setFloat32(24, h8, true);
        globalsView.setUint32(28, this.numParticles, true);
        this.queue.writeBuffer(this.sphCalcUniformBuffer, 0, globalsBuffer);

        const boundsGlobalUniformBuffer = new Float32Array([
            bounds.minX,
            bounds.minY,
            bounds.maxX,
            bounds.maxY,
            this.gridWidth,
            this.gridHeight,
            0,0
        ])
        this.queue.writeBuffer(this.boundsUniformBuffer, 0, boundsGlobalUniformBuffer);

        const sphCommandEncoder = this.device.createCommandEncoder();
        const passEncoder = sphCommandEncoder.beginComputePass();

        passEncoder.setPipeline(this.computePipeline)
        passEncoder.setBindGroup(0, this.sphCalcBindGroup);

        const workgroupCount = Math.ceil(this.numParticles / 64);
        passEncoder.dispatchWorkgroups(workgroupCount);
        passEncoder.end();

        sphCommandEncoder.copyBufferToBuffer(
            this.positionsBuffer,
            0,
            this.positionsReadbackBuffer,
            0,
            this.positionsBufferSize
        );

        sphCommandEncoder.copyBufferToBuffer(
            this.densitiesBuffer,
            0,
            this.densitiesReadbackBuffer,
            0,
            Float32Array.BYTES_PER_ELEMENT * this.numParticles
        );

        this.queue.submit([sphCommandEncoder.finish()]);
        await this.queue.onSubmittedWorkDone();
        
        await this.readBufferIntoArray(this.positionsReadbackBuffer, this.positionsBufferSize, this.positions);
        await this.readBufferIntoArray(this.densitiesReadbackBuffer, Float32Array.BYTES_PER_ELEMENT * this.numParticles, this.densities);
        
        this.#particles.forEach(particle => {
            particle.position.x = this.positions[particle.index * 2];
            particle.position.y = this.positions[(particle.index * 2) + 1];
        });
        if(this.forceDirection){
            this.applyOutwardForce(bounds);
        }else{
            this.applyInwardForce(bounds);
        }
    }

    async readBufferIntoArray(buffer, byteLength, targetArray) {
        await buffer.mapAsync(GPUMapMode.READ, 0, byteLength);
        const mapped = buffer.getMappedRange(0, byteLength);
        targetArray.set(new Float32Array(mapped));
        buffer.unmap();
    }

    async update(now, bounds) {
        if (this.lastTime === 0) {
            this.lastTime = now;
            // Initialize grid on first frame
            this.rebuildGrid(bounds);
            return;
        }

        let frameDelta = (now - this.lastTime) / 1000;
        this.lastTime = now;
        frameDelta = Math.min(frameDelta, this.MAX_ACCUM);

        this.accumulator += frameDelta;

        let steps = 0;
        while (this.accumulator >= this.FIXED_DT && steps < this.MAX_STEPS) {
            // await this.applyWgslPhysics(this.FIXED_DT, bounds);
            this.stepPhysics(this.FIXED_DT,bounds)
            this.accumulator -= this.FIXED_DT;
            steps++;
        }

        this.setColorBasedOnDensity();
    }

    getParticles() {
        return this.#particles;
    }

    applyOutwardForce(bounds, radius = 5.0, strength = 125) {
        if (this.mouseX === undefined || this.mouseY === undefined) return;

        const rMax2 = radius * radius;

        const cellX = Math.floor((this.mouseX - bounds.minX) / this.H);
        const cellY = Math.floor((this.mouseY - bounds.minY) / this.H);
        if (cellX < 0 || cellX >= this.gridWidth || cellY < 0 || cellY >= this.gridHeight) return;

        const neighborCells = this.getNeighborCellsInRadius(cellX, cellY, radius);

        for (let c = 0; c < neighborCells.length; c++) {
            const cell = neighborCells[c];
            for (let n = 0; n < cell.length; n++) {
                const pIdx = cell[n];
                const pX = this.positions[pIdx * 2];
                const pY = this.positions[pIdx * 2 + 1];

                const dx = pX - this.mouseX;
                const dy = pY - this.mouseY;
                const r2 = dx * dx + dy * dy;

                if (r2 === 0 || r2 > rMax2) continue;

                const r = Math.sqrt(r2);
                const invR = 1 / r;
                const falloff = 1 - (r / radius);
                const scaled = strength * falloff * falloff;

                this.accelerations[pIdx * 2]     += dx * invR * scaled;
                this.accelerations[pIdx * 2 + 1] += dy * invR * scaled;
            }
        }
    }

    applyInwardForce(bounds, radius = 10.0, strength = 125) {
        if (this.mouseX === undefined || this.mouseY === undefined) return;

        const rMax2 = radius * radius;

        const cellX = Math.floor((this.mouseX - bounds.minX) / this.H);
        const cellY = Math.floor((this.mouseY - bounds.minY) / this.H);
        if (cellX < 0 || cellX >= this.gridWidth || cellY < 0 || cellY >= this.gridHeight) return;

        const neighborCells = this.getNeighborCellsInRadius(cellX, cellY, radius);

        for (let c = 0; c < neighborCells.length; c++) {
            const cell = neighborCells[c];
            for (let n = 0; n < cell.length; n++) {
                const pIdx = cell[n];
                const pX = this.positions[pIdx * 2];
                const pY = this.positions[pIdx * 2 + 1];

                const dx = pX - this.mouseX;
                const dy = pY - this.mouseY;
                const r2 = dx * dx + dy * dy;

                if (r2 === 0 || r2 > rMax2) continue;

                const r = Math.sqrt(r2);
                const invR = 1 / r;
                const falloff = 1 - (r / radius);
                const scaled = strength * falloff * falloff;

                this.accelerations[pIdx * 2]     -= dx * invR * scaled;
                this.accelerations[pIdx * 2 + 1] -= dy * invR * scaled;
            }
        }
    }

    addParticles(num, scene) {
        const oldCount = this.numParticles;
        const newCount = oldCount + num;

        const newPositions = new Float32Array(newCount * 2);
        const newVelocities = new Float32Array(newCount * 2);
        const newAccelerations = new Float32Array(newCount * 2);
        const newPressures = new Float32Array(newCount);
        const newDensities = new Float32Array(newCount);

        newPositions.set(this.positions);
        newVelocities.set(this.velocities);
        newAccelerations.set(this.accelerations);
        newPressures.set(this.pressures);
        newDensities.set(this.densities);

        this.positions = newPositions;
        this.velocities = newVelocities;
        this.accelerations = newAccelerations;
        this.pressures = newPressures;
        this.densities = newDensities;

        for (let i = 0; i < num; i++) {
            const index = oldCount + i;
            
            const geometry = new THREE.CircleGeometry(0.07, 10);
            const material = new THREE.MeshBasicMaterial({ color: 0x4fd9f7 });
            let newParticle = new particle(geometry, material);
            
            newParticle.index = index;
            newParticle.position.x = this.mouseX || 0;
            newParticle.position.y = this.mouseY || 0;
            
            this.positions[index * 2] = newParticle.position.x;
            this.positions[index * 2 + 1] = newParticle.position.y;
            this.velocities[index * 2] = 0;
            this.velocities[index * 2 + 1] = 0;
            this.accelerations[index * 2] = 0;
            this.accelerations[index * 2 + 1] = this.GRAVITY_Y;
            this.densities[index] = 0;
            this.pressures[index] = 0;
            
            this.#particles.push(newParticle);
            scene.add(newParticle);
        }
        
        this.numParticles = this.#particles.length;
    }

    removeParticles(num, scene) {
        const actualRemove = Math.min(num, this.numParticles);
        const newCount = this.numParticles - actualRemove;

        for (let i = 0; i < actualRemove; i++) {
            const removedParticle = this.#particles.pop();
            scene.remove(removedParticle);
        }

        this.positions = this.positions.slice(0, newCount * 2);
        this.velocities = this.velocities.slice(0, newCount * 2);
        this.accelerations = this.accelerations.slice(0, newCount * 2);
        this.pressures = this.pressures.slice(0, newCount);
        this.densities = this.densities.slice(0, newCount);

        this.numParticles = this.#particles.length;
    }

    setMousePosWorldCord(x,y){
        this.mouseX = x;
        this.mouseY = y;
    }

    setRestDensity(rest){
        this.REST_DENSITY = parseFloat(rest);
    }

    setGravity(value){
        this.GRAVITY_Y = -parseFloat(value);
    }

    setGasConst(value){
        this.GAS_CONST = parseFloat(value);
    }

    setViscocity(value){
        this.VISCOSITY_COEFF = parseFloat(value);
    }
}

export default Engine;