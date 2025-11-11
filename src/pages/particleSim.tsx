import { useEffect, useRef } from "react";
import SimMenu from "../components/simMenu";

function ParticleSim() {
  const canvasHostRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    let dispose: (() => void) | undefined;
    let mounted = true;

    (async () => {
      try {
        const module = await import("../particle-sim-js/main.js");
        if (!mounted || !canvasHostRef.current) {
          return;
        }
        dispose = await module.startParticleSim(canvasHostRef.current);
      } catch (error) {
        console.error("Failed to start particle simulation", error);
      }
    })();

    return () => {
      mounted = false;
      dispose?.();
      dispose = undefined;
    };
  }, []);

  return (
    <div style={{ position: "relative", minHeight: "100vh" }}>
      <div
        ref={canvasHostRef}
        style={{ position: "fixed", inset: 0, zIndex: 0 }}
      />
      <SimMenu />
    </div>
  );
}

export default ParticleSim;