import * as THREE from 'three/webgpu';

class particle extends THREE.Mesh {
    index;
    mass = 1;
    constructor(geometry, material) {
        super(geometry,material);
    }
}

class Engine {
  instancedMesh = null;
  #dummy = new THREE.Object3D();
  #color = new THREE.Color();
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
  H = 0.3; // also grid cell size
  GAS_CONST = 80;
  GRAVITY_Y = -9;
  DENSITY_EPS = 1e-6; // to avoid divide-by-zero
  PARTICLE_RADIUS = 0.03;

  RESTITUTION = 0.5;
  WALL_FRICTION = 0.1;
  REST_DENSITY = 4;
  VISCOSITY_COEFF = 80;

  // Fixed-step physics update for stability
  FIXED_DT = 1 / 100; // seconds per physics step
  MAX_ACCUM = 0.05; // clamp big frame gaps (s)
  MAX_STEPS = 8; // prevent spiral of death

  lastTime = 0;
  accumulator = 0;


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
    const gradMagnitude = (constant * term) / dist;

    return {
      gradWx: gradMagnitude * dx,
      gradWy: gradMagnitude * dy,
    };
  }



  getCellCoords(particleIndex, bounds) {
    const px = this.positions[particleIndex * 2];
    const py = this.positions[particleIndex * 2 + 1];
    const cellX = Math.floor((px - bounds.minX) / this.H);
    const cellY = Math.floor((py - bounds.minY) / this.H);
    return { cellX, cellY };
  }

  createParticles(numParticles, minX, maxX, minY, maxY, scene) {
    const maxCapacity = 10000;
    this.numParticles = numParticles;

    this.positions = new Float32Array(numParticles * 2);
    this.velocities = new Float32Array(numParticles * 2);
    this.accelerations = new Float32Array(numParticles * 2);

    this.pressures = new Float32Array(numParticles);
    this.viscocities = new Float32Array(numParticles);
    this.densities = new Float32Array(numParticles);

    // create grid
    this.gridWidth = Math.ceil((maxX - minX) / this.H);
    this.gridHeight = Math.ceil((maxY - minY) / this.H);
    const totalCells = this.gridWidth * this.gridHeight;

    for (let cellIndex = 0; cellIndex < totalCells; cellIndex++) {
      this.cells.push([]);
    }

    // even distribution - grid layout
    const gridCols = Math.ceil(Math.sqrt(numParticles));
    const gridRows = Math.ceil(numParticles / gridCols);
    const spacingX = (maxX - minX) / (gridCols + 1);
    const spacingY = (maxY - minY) / (gridRows + 1);

    const geometry = new THREE.CircleGeometry(0.07, 10);
    const material = new THREE.MeshBasicMaterial({ color: 0xffffff });
    this.instancedMesh = new THREE.InstancedMesh(geometry, material, maxCapacity);
    this.instancedMesh.count = numParticles;
    this.instancedMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);

    for (let index = 0; index < numParticles; index++) {
      const col = index % gridCols;
      const row = Math.floor(index / gridCols);

      const px = minX + (col + 1) * spacingX;
      const py = minY + (row + 1) * spacingY;

      this.positions[index * 2] = px;
      this.positions[index * 2 + 1] = py;

      const cellY = Math.floor((py - minY) / this.H);
      const cellX = Math.floor((px - minX) / this.H);
      this.cells[cellY * this.gridWidth + cellX].push(index);

      this.#dummy.position.set(px, py, 0);
      this.#dummy.updateMatrix();
      this.instancedMesh.setMatrixAt(index, this.#dummy.matrix);
      this.instancedMesh.setColorAt(index, this.#color.set(0x4fd9f7));
    }

    this.instancedMesh.instanceMatrix.needsUpdate = true;
    if (this.instancedMesh.instanceColor) {
      this.instancedMesh.instanceColor.setUsage(THREE.DynamicDrawUsage);
      this.instancedMesh.instanceColor.needsUpdate = true;
    }
    scene.add(this.instancedMesh);
    return this.instancedMesh;
  }

  applyPhysics(dt) {
    // Semi-implicit Euler with velocity clamping
    const MAX_SPEED = 25.0;
    const MAX_SPEED_SQ = MAX_SPEED * MAX_SPEED;

    for (let idx = 0; idx < this.numParticles; idx++) {
      let vx = this.velocities[idx * 2] + this.accelerations[idx * 2] * dt;
      let vy = this.velocities[idx * 2 + 1] + this.accelerations[idx * 2 + 1] * dt;

      const speedSq = vx * vx + vy * vy;
      if (speedSq > MAX_SPEED_SQ) {
        const factor = MAX_SPEED / Math.sqrt(speedSq);
        vx *= factor;
        vy *= factor;
      }

      this.velocities[idx * 2] = vx;
      this.velocities[idx * 2 + 1] = vy;

      this.positions[idx * 2] += vx * dt;
      this.positions[idx * 2 + 1] += vy * dt;
    }
  }

  calculateDensity(bounds) {
    const h2 = this.H * this.H;
    const EPSILON = 1e-3;
    const EPSILON_SQ = EPSILON * EPSILON;

    for (let targetElem = 0; targetElem < this.numParticles; targetElem++) {
      let tx = this.positions[targetElem * 2];
      let ty = this.positions[targetElem * 2 + 1];
      const cellX = Math.floor((tx - bounds.minX) / this.H);
      const cellY = Math.floor((ty - bounds.minY) / this.H);

      const minCX = Math.max(0, cellX - 1);
      const maxCX = Math.min(this.gridWidth - 1, cellX + 1);
      const minCY = Math.max(0, cellY - 1);
      const maxCY = Math.min(this.gridHeight - 1, cellY + 1);

      let rho = 0;
      for (let cy = minCY; cy <= maxCY; cy++) {
        const rowOffset = cy * this.gridWidth;
        for (let cx = minCX; cx <= maxCX; cx++) {
          const cell = this.cells[rowOffset + cx];
          for (let k = 0; k < cell.length; k++) {
            const outerElem = cell[k];
            if (outerElem === targetElem) continue;

            const ox = this.positions[outerElem * 2];
            const oy = this.positions[outerElem * 2 + 1];
            let dx = ox - tx;
            let dy = oy - ty;
            let r2 = dx * dx + dy * dy;

            // Separate overlapping particles
            if (r2 < EPSILON_SQ) {
              const angle = Math.random() * Math.PI * 2;
              const offX = Math.cos(angle) * EPSILON * 0.5;
              const offY = Math.sin(angle) * EPSILON * 0.5;

              tx += offX;
              ty += offY;
              this.positions[targetElem * 2] = tx;
              this.positions[targetElem * 2 + 1] = ty;
              this.positions[outerElem * 2] -= offX;
              this.positions[outerElem * 2 + 1] -= offY;

              dx = (ox - offX) - tx;
              dy = (oy - offY) - ty;
              r2 = dx * dx + dy * dy;
            }

            if (r2 <= h2) {
              rho += this.poly6Kernel(r2, this.H);
            }
          }
        }
      }
      this.densities[targetElem] = rho;
    }
  }

  calculatePressure() {
    for (let i = 0; i < this.numParticles; i++) {
      const pressure =
        this.GAS_CONST * (this.densities[i] - this.REST_DENSITY);
      this.pressures[i] = Math.max(0, pressure);
    }
  }

  calculatePressureForce(bounds) {
    const h2 = this.H * this.H;

    for (let targetElem = 0; targetElem < this.numParticles; targetElem++) {
      const ix = this.positions[targetElem * 2];
      const iy = this.positions[targetElem * 2 + 1];
      const cellX = Math.floor((ix - bounds.minX) / this.H);
      const cellY = Math.floor((iy - bounds.minY) / this.H);

      const minCX = Math.max(0, cellX - 1);
      const maxCX = Math.min(this.gridWidth - 1, cellX + 1);
      const minCY = Math.max(0, cellY - 1);
      const maxCY = Math.min(this.gridHeight - 1, cellY + 1);

      let fx = 0,
        fy = 0;

      for (let cy = minCY; cy <= maxCY; cy++) {
        const rowOffset = cy * this.gridWidth;
        for (let cx = minCX; cx <= maxCX; cx++) {
          const cell = this.cells[rowOffset + cx];
          for (let k = 0; k < cell.length; k++) {
            const outerElem = cell[k];
            if (outerElem === targetElem) {
              continue;
            }

            const jx = this.positions[outerElem * 2];
            const jy = this.positions[outerElem * 2 + 1];

            const dx = ix - jx;
            const dy = iy - jy;
            const r2 = dx * dx + dy * dy;

            if (r2 > 0 && r2 <= h2) {
              const { gradWx, gradWy } = this.spikyKernel(dx, dy, r2, this.H);

              const density_i = Math.max(
                this.densities[targetElem],
                this.DENSITY_EPS,
              );
              const density_j = Math.max(
                this.densities[outerElem],
                this.DENSITY_EPS,
              );
              const denom_i = density_i * density_i;
              const denom_j = density_j * density_j;

              const pressure_i = this.pressures[targetElem];
              const pressure_j = this.pressures[outerElem];

              const coeff = pressure_i / denom_i + pressure_j / denom_j;

              fx -= coeff * gradWx;
              fy -= coeff * gradWy;
            }
          }
        }
      }

      this.accelerations[targetElem * 2] += fx;
      this.accelerations[targetElem * 2 + 1] += fy;
    }
  }

  applyViscosity(dt, bounds) {
    const h2 = this.H * this.H;

    for (let targetElem = 0; targetElem < this.numParticles; targetElem++) {
      const ix = this.positions[targetElem * 2];
      const iy = this.positions[targetElem * 2 + 1];
      const cellX = Math.floor((ix - bounds.minX) / this.H);
      const cellY = Math.floor((iy - bounds.minY) / this.H);

      const minCX = Math.max(0, cellX - 1);
      const maxCX = Math.min(this.gridWidth - 1, cellX + 1);
      const minCY = Math.max(0, cellY - 1);
      const maxCY = Math.min(this.gridHeight - 1, cellY + 1);

      let corrX = 0,
        corrY = 0;

      for (let cy = minCY; cy <= maxCY; cy++) {
        const rowOffset = cy * this.gridWidth;
        for (let cx = minCX; cx <= maxCX; cx++) {
          const cell = this.cells[rowOffset + cx];
          for (let k = 0; k < cell.length; k++) {
            const outerElem = cell[k];
            if (outerElem === targetElem) {
              continue;
            }

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
              const densityJ = Math.max(
                this.densities[outerElem],
                this.DENSITY_EPS,
              );
              corrX += (jvx - ivx) * (1 / densityJ) * w;
              corrY += (jvy - ivy) * (1 / densityJ) * w;
            }
          }
        }
      }

      this.velocities[targetElem * 2] += this.VISCOSITY_COEFF * corrX * dt;
      this.velocities[targetElem * 2 + 1] += this.VISCOSITY_COEFF * corrY * dt;
    }
  }

    setColorBasedOnDensity() {
        let minD = Infinity;
        let maxD = -Infinity;
        
        for (let i = 0; i < this.numParticles; i++) {
            if (this.densities[i] < minD) minD = this.densities[i];
            if (this.densities[i] > maxD) maxD = this.densities[i];
        }
        
        const range = Math.max(maxD - minD, 1e-8);

        for (let i = 0; i < this.numParticles; i++) {
            const t = (this.densities[i] - minD) / range;
            const red = t * 1.7;
            const green = 0;
            const blue = 1 - t;
            this.#color.setRGB(red, green, blue);
            this.instancedMesh.setColorAt(i, this.#color);
        }
        if (this.instancedMesh.instanceColor) {
            this.instancedMesh.instanceColor.needsUpdate = true;
        }
    }

  handleCollisions(bounds) {
    const minX = bounds.minX + this.PARTICLE_RADIUS;
    const maxX = bounds.maxX - this.PARTICLE_RADIUS;
    const minY = bounds.minY + this.PARTICLE_RADIUS;
    const maxY = bounds.maxY - this.PARTICLE_RADIUS;

    for (let idx = 0; idx < this.numParticles; idx++) {
      let px = this.positions[idx * 2];
      let py = this.positions[idx * 2 + 1];
      let vx = this.velocities[idx * 2];
      let vy = this.velocities[idx * 2 + 1];

      // Left/right walls
      if (px < minX) {
        px = minX;
        if (vx < 0) vx = -vx * this.RESTITUTION;
        vy *= this.WALL_FRICTION;
      } else if (px > maxX) {
        px = maxX;
        if (vx > 0) vx = -vx * this.RESTITUTION;
        vy *= this.WALL_FRICTION;
      }

      // Bottom/top walls
      if (py < minY) {
        py = minY;
        if (vy < 0) vy = -vy * this.RESTITUTION;
        vx *= this.WALL_FRICTION;
      } else if (py > maxY) {
        py = maxY;
        if (vy > 0) vy = -vy * this.RESTITUTION;
        vx *= this.WALL_FRICTION;
      }

      // Write back to arrays
      this.positions[idx * 2] = px;
      this.positions[idx * 2 + 1] = py;
      this.velocities[idx * 2] = vx;
      this.velocities[idx * 2 + 1] = vy;
    }
  }

  rebuildGrid(bounds) {
    const nextWidth = Math.ceil((bounds.maxX - bounds.minX) / this.H);
    const nextHeight = Math.ceil((bounds.maxY - bounds.minY) / this.H);
    const totalCells = nextWidth * nextHeight;

    if (
      this.gridWidth === nextWidth &&
      this.gridHeight === nextHeight &&
      this.cells.length === totalCells
    ) {
      return;
    }

    this.gridWidth = nextWidth;
    this.gridHeight = nextHeight;

    this.cells = [];
    for (let cellIndex = 0; cellIndex < totalCells; cellIndex++) {
      this.cells.push([]);
    }
  }

  stepPhysics(dt, bounds) {
    // Reset accelerations
    for (let i = 0; i < this.accelerations.length / 2; i++) {
      this.accelerations[i * 2] = 0;
      this.accelerations[i * 2 + 1] = this.GRAVITY_Y;
    }

    this.rebuildGrid(bounds);
    for (let c = 0; c < this.cells.length; c++) {
      this.cells[c].length = 0;
    }
    for (let i = 0; i < this.numParticles; i++) {
      const px = this.positions[i * 2];
      const py = this.positions[i * 2 + 1];
      const cellX = Math.floor((px - bounds.minX) / this.H);
      const cellY = Math.floor((py - bounds.minY) / this.H);
      if (
        cellX >= 0 &&
        cellX < this.gridWidth &&
        cellY >= 0 &&
        cellY < this.gridHeight
      ) {
        this.cells[cellY * this.gridWidth + cellX].push(i);
      }
    }

    // SPH pipeline
    this.calculateDensity(bounds);
    this.calculatePressure();
    this.calculatePressureForce(bounds);
    this.applyViscosity(dt, bounds);
    if (this.forceDirection) {
      this.applyOutwardForce(bounds);
    } else {
      this.applyInwardForce(bounds);
    }

    this.applyPhysics(dt);
  }

  update(now, bounds) {
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
      this.stepPhysics(this.FIXED_DT, bounds);
      this.handleCollisions(bounds);
      this.accumulator -= this.FIXED_DT;
      steps++;
    }

    this.setColorBasedOnDensity();

    if (this.instancedMesh) {
      for (let i = 0; i < this.numParticles; i++) {
        this.#dummy.position.set(
          this.positions[i * 2],
          this.positions[i * 2 + 1],
          0,
        );
        this.#dummy.updateMatrix();
        this.instancedMesh.setMatrixAt(i, this.#dummy.matrix);
      }
      this.instancedMesh.instanceMatrix.needsUpdate = true;
    }
  }

  getParticles() {
    return [];
  }

  applyOutwardForce(bounds, radius = 5.0, strength = 125) {
    if (this.mouseX === undefined || this.mouseY === undefined) return;

    const rMax2 = radius * radius;

    const cellX = Math.floor((this.mouseX - bounds.minX) / this.H);
    const cellY = Math.floor((this.mouseY - bounds.minY) / this.H);
    if (
      cellX < 0 ||
      cellX >= this.gridWidth ||
      cellY < 0 ||
      cellY >= this.gridHeight
    )
      return;

    const cellR = Math.ceil(radius / this.H);
    const minCX = Math.max(0, cellX - cellR);
    const maxCX = Math.min(this.gridWidth - 1, cellX + cellR);
    const minCY = Math.max(0, cellY - cellR);
    const maxCY = Math.min(this.gridHeight - 1, cellY + cellR);

    for (let y = minCY; y <= maxCY; y++) {
      const rowOffset = y * this.gridWidth;
      for (let x = minCX; x <= maxCX; x++) {
        const cell = this.cells[rowOffset + x];
        for (let n = 0; n < cell.length; n++) {
          const pIdx = cell[n];
          const pX = this.positions[pIdx * 2];
          const pY = this.positions[pIdx * 2 + 1];

          const dx = pX - this.mouseX;
          const dy = pY - this.mouseY;
          const r2 = dx * dx + dy * dy;

          if (r2 === 0 || r2 > rMax2) continue;

          const r = Math.sqrt(r2);
          const invR = 1 / Math.max(r, 0.1);
          const falloff = 1 - r / radius;
          const scaled = strength * falloff * falloff;

          this.accelerations[pIdx * 2] += dx * invR * scaled;
          this.accelerations[pIdx * 2 + 1] += dy * invR * scaled;
        }
      }
    }
  }

  applyInwardForce(bounds, radius = 10.0, strength = 125) {
    if (this.mouseX === undefined || this.mouseY === undefined) return;

    const rMax2 = radius * radius;

    const cellX = Math.floor((this.mouseX - bounds.minX) / this.H);
    const cellY = Math.floor((this.mouseY - bounds.minY) / this.H);
    if (
      cellX < 0 ||
      cellX >= this.gridWidth ||
      cellY < 0 ||
      cellY >= this.gridHeight
    )
      return;

    const cellR = Math.ceil(radius / this.H);
    const minCX = Math.max(0, cellX - cellR);
    const maxCX = Math.min(this.gridWidth - 1, cellX + cellR);
    const minCY = Math.max(0, cellY - cellR);
    const maxCY = Math.min(this.gridHeight - 1, cellY + cellR);

    for (let y = minCY; y <= maxCY; y++) {
      const rowOffset = y * this.gridWidth;
      for (let x = minCX; x <= maxCX; x++) {
        const cell = this.cells[rowOffset + x];
        for (let n = 0; n < cell.length; n++) {
          const pIdx = cell[n];
          const pX = this.positions[pIdx * 2];
          const pY = this.positions[pIdx * 2 + 1];

          const dx = pX - this.mouseX;
          const dy = pY - this.mouseY;
          const r2 = dx * dx + dy * dy;

          if (r2 === 0 || r2 > rMax2) continue;

          const r = Math.sqrt(r2);
          const invR = 1 / Math.max(r, 0.1);
          const falloff = 1 - r / radius;
          const scaled = strength * falloff * falloff;

          this.accelerations[pIdx * 2] -= dx * invR * scaled;
          this.accelerations[pIdx * 2 + 1] -= dy * invR * scaled;
        }
      }
    }
  }

  addParticles(num, scene) {
    const oldCount = this.numParticles;
    const newCount = Math.min(10000, oldCount + num);
    const actualAdd = newCount - oldCount;
    if (actualAdd <= 0) return;

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

    for (let i = 0; i < actualAdd; i++) {
      const index = oldCount + i;
      this.positions[index * 2] = (this.mouseX || 0) + (Math.random() - 0.5) * 0.05;
      this.positions[index * 2 + 1] = (this.mouseY || 0) + (Math.random() - 0.5) * 0.05;
      this.velocities[index * 2] = 0;
      this.velocities[index * 2 + 1] = 0;
      this.accelerations[index * 2] = 0;
      this.accelerations[index * 2 + 1] = this.GRAVITY_Y;
      this.densities[index] = 0;
      this.pressures[index] = 0;
    }

    this.numParticles = newCount;
    if (this.instancedMesh) {
      this.instancedMesh.count = newCount;
    }
  }

  removeParticles(num, scene) {
    const actualRemove = Math.min(num, this.numParticles);
    const newCount = this.numParticles - actualRemove;

    this.positions = this.positions.slice(0, newCount * 2);
    this.velocities = this.velocities.slice(0, newCount * 2);
    this.accelerations = this.accelerations.slice(0, newCount * 2);
    this.pressures = this.pressures.slice(0, newCount);
    this.densities = this.densities.slice(0, newCount);

    this.numParticles = newCount;
    if (this.instancedMesh) {
      this.instancedMesh.count = newCount;
    }
  }

  dispose(scene) {
    if (this.instancedMesh) {
      if (scene) {
        scene.remove(this.instancedMesh);
      }
      this.instancedMesh.geometry.dispose();
      this.instancedMesh.material.dispose();
      this.instancedMesh.dispose();
      this.instancedMesh = null;
    }
  }

  setMousePosWorldCord(x, y) {
    this.mouseX = x;
    this.mouseY = y;
  }

  setRestDensity(rest) {
    this.REST_DENSITY = parseFloat(rest);
  }

  setGravity(value) {
    this.GRAVITY_Y = -parseFloat(value);
  }

  setGasConst(value) {
    this.GAS_CONST = parseFloat(value);
  }

  setViscocity(value) {
    this.VISCOSITY_COEFF = parseFloat(value);
  }
}

export default Engine;