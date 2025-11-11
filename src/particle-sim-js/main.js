import * as THREE from 'three/webgpu';
import Engine from './engine.js';

const raycaster = new THREE.Raycaster();
const planeZ0 = new THREE.Plane(new THREE.Vector3(0, 0, 1), 0);
const tmpPt = new THREE.Vector3();

let renderer = null;
let scene = null;
let camera = null;
let engine = null;
let bounds = null;
let canvasHost = null;
let pointerDown = false;
let animationFrameId = 0;

let resizeHandler = null;
let pointerDownHandler = null;
let pointerMoveHandler = null;
let pointerUpHandler = null;
let pointerLeaveHandler = null;
let pointerCancelHandler = null;

function calculateViewportBounds(currentCamera) {
  const fov = currentCamera.fov * (Math.PI / 180);
  const height = 2 * Math.tan(fov / 2) * Math.abs(currentCamera.position.z);
  const width = height * currentCamera.aspect;
  return { minX: -width / 2, maxX: width / 2, minY: -height / 2, maxY: height / 2 };
}

function setupScene(hostElement) {
  canvasHost = hostElement;
  scene = new THREE.Scene();

  const width = hostElement.clientWidth || window.innerWidth;
  const height = hostElement.clientHeight || window.innerHeight;

  camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
  renderer = new THREE.WebGPURenderer();
  renderer.setSize(width, height);
  renderer.domElement.style.width = "100%";
  renderer.domElement.style.height = "100%";
  renderer.domElement.style.display = "block";
  hostElement.appendChild(renderer.domElement);

  scene.background = new THREE.Color("black");
  camera.position.z = 10;

  resizeHandler = () => {
    if (!renderer || !camera) {
      return;
    }
    const nextWidth = hostElement.clientWidth || window.innerWidth;
    const nextHeight = hostElement.clientHeight || window.innerHeight;
    camera.aspect = nextWidth / nextHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(nextWidth, nextHeight);
    bounds = calculateViewportBounds(camera);
  };

  pointerDownHandler = (event) => {
    if (!engine || event.button !== 0) {
      return;
    }
    if (event.target instanceof HTMLElement && event.target.closest(".settings-menu")) {
      return;
    }
    pointerDown = true;
    const wp = screenPosToWorldPos(event.clientX, event.clientY);
    if (wp) {
      engine.setMousePosWorldCord(wp.x, wp.y);
    }
  };

  pointerMoveHandler = (event) => {
    if (!engine || !pointerDown) {
      return;
    }
    const wp = screenPosToWorldPos(event.clientX, event.clientY);
    if (wp) {
      engine.setMousePosWorldCord(wp.x, wp.y);
    }
  };

  const clearPointer = () => {
    if (!engine) {
      return;
    }
    pointerDown = false;
    engine.setMousePosWorldCord(Number.NaN, Number.NaN);
  };

  pointerUpHandler = clearPointer;
  pointerLeaveHandler = clearPointer;
  pointerCancelHandler = clearPointer;

  window.addEventListener("resize", resizeHandler);
  window.addEventListener("pointerdown", pointerDownHandler);
  window.addEventListener("pointermove", pointerMoveHandler);
  window.addEventListener("pointerup", pointerUpHandler);
  window.addEventListener("pointerleave", pointerLeaveHandler);
  window.addEventListener("pointercancel", pointerCancelHandler);
}

function screenPosToWorldPos(clientX, clientY) {
  if (!renderer || !camera) {
    return null;
  }

  const rect = renderer.domElement.getBoundingClientRect();
  const ndcX = ((clientX - rect.left) / rect.width) * 2 - 1;
  const ndcY = -((clientY - rect.top) / rect.height) * 2 + 1;

  raycaster.setFromCamera({ x: ndcX, y: ndcY }, camera);
  if (raycaster.ray.intersectPlane(planeZ0, tmpPt)) {
    return tmpPt.clone();
  }

  return null;
}

function attachSliderHandlers() {
  const restDensitySlider = document.getElementById("density");
  const gravitySlider = document.getElementById("gravity");
  const gasConstSlider = document.getElementById("gas-const");
  const viscocitySlider = document.getElementById("viscocity");
  const particleSlider = document.getElementById("particles");
  const mouseBehaviour = document.getElementById("mouse-behaviour");

  if (
    !engine ||
    !(restDensitySlider instanceof HTMLInputElement) ||
    !(gravitySlider instanceof HTMLInputElement) ||
    !(gasConstSlider instanceof HTMLInputElement) ||
    !(viscocitySlider instanceof HTMLInputElement) ||
    !(particleSlider instanceof HTMLInputElement)
  ) {
    throw new Error("Simulation controls are missing from the DOM.");
  }

  const disposers = [];

  const mouseBehaviourHandler = () => {
    engine.forceDirection = !engine.forceDirection;
  };
  mouseBehaviour.addEventListener("change", mouseBehaviourHandler);
  disposers.push(() => mouseBehaviour.removeEventListener("change", mouseBehaviourHandler));

  const restDensityHandler = (event) => engine.setRestDensity(event.target.value);
  restDensitySlider.addEventListener("input", restDensityHandler);
  disposers.push(() => restDensitySlider.removeEventListener("input", restDensityHandler));

  const gravityHandler = (event) => engine.setGravity(event.target.value);
  gravitySlider.addEventListener("input", gravityHandler);
  disposers.push(() => gravitySlider.removeEventListener("input", gravityHandler));

  const gasConstHandler = (event) => engine.setGasConst(event.target.value);
  gasConstSlider.addEventListener("input", gasConstHandler);
  disposers.push(() => gasConstSlider.removeEventListener("input", gasConstHandler));

  const viscocityHandler = (event) => engine.setViscocity(event.target.value);
  viscocitySlider.addEventListener("input", viscocityHandler);
  disposers.push(() => viscocitySlider.removeEventListener("input", viscocityHandler));

  const particleHandler = (event) => {
    const desiredCount = Number(event.target.value);
    const difference = desiredCount - engine.numParticles;
    if (difference > 0) {
      engine.addParticles(difference, scene);
    } else if (difference < 0) {
      engine.removeParticles(-difference, scene);
    }
  };
  particleSlider.addEventListener("input", particleHandler);
  disposers.push(() => particleSlider.removeEventListener("input", particleHandler));

  return () => {
    disposers.splice(0).forEach((dispose) => dispose());
  };
}

function animate(now) {
  if (!engine || !renderer || !camera || !bounds) {
    return;
  }
  animationFrameId = requestAnimationFrame(animate);
  engine.update(now, bounds);
  renderer.renderAsync(scene, camera);
}

export async function startParticleSim(hostElement) {
  engine = new Engine();
  await engine.init();

  const targetHost = hostElement || document.body;
  setupScene(targetHost);
  bounds = calculateViewportBounds(camera);

  const sliderCleanup = attachSliderHandlers();

  const particleCount = 3000;
  engine.createParticles(
    particleCount,
    bounds.minX,
    bounds.maxX,
    bounds.minY,
    bounds.maxY,
    scene
  );

  animationFrameId = requestAnimationFrame(animate);

  return () => {
    sliderCleanup();
    cancelAnimationFrame(animationFrameId);
    animationFrameId = 0;

    if (resizeHandler) {
      window.removeEventListener("resize", resizeHandler);
      resizeHandler = null;
    }
    if (pointerDownHandler) {
      window.removeEventListener("pointerdown", pointerDownHandler);
      pointerDownHandler = null;
    }
    if (pointerMoveHandler) {
      window.removeEventListener("pointermove", pointerMoveHandler);
      pointerMoveHandler = null;
    }
    if (pointerUpHandler) {
      window.removeEventListener("pointerup", pointerUpHandler);
      pointerUpHandler = null;
    }
    if (pointerLeaveHandler) {
      window.removeEventListener("pointerleave", pointerLeaveHandler);
      pointerLeaveHandler = null;
    }
    if (pointerCancelHandler) {
      window.removeEventListener("pointercancel", pointerCancelHandler);
      pointerCancelHandler = null;
    }

    pointerDown = false;

    if (engine) {
      engine.setMousePosWorldCord(Number.NaN, Number.NaN);
    }

    if (renderer) {
      renderer.dispose();
      if (renderer.domElement && renderer.domElement.parentElement === targetHost) {
        targetHost.removeChild(renderer.domElement);
      }
    }

    renderer = null;
    scene = null;
    camera = null;
    bounds = null;
    canvasHost = null;
    engine = null;
  };
}
