declare module "*/particle-sim-js/main.js" {
  export function startParticleSim(
    hostElement: HTMLElement
  ): Promise<() => void>;
}