import { useState } from "react";
import "../css/simMenu.css";

function SimMenu() {
  const [restDensity, setRestDensity] = useState(4);
  const [gasConst, setGasConst] = useState(80);
  const [gravity, setGravity] = useState(9);
  const [viscosity, setViscosity] = useState(80);
  const [particleCount, setParticleCount] = useState(3000);
  const [mouseEffect, setMouseEffect] = useState(true);
  const [isOpen, setIsOpen] = useState(true);

  return (
    <div
      className={`settings-wrapper ${
        isOpen ? "settings-wrapper--open" : "settings-wrapper--closed"
      }`}
    >
      <button
        type="button"
        className="settings-tab"
        aria-expanded={isOpen}
        aria-controls="sim-menu-panel"
        onClick={() => setIsOpen((prev) => !prev)}
      >
        {isOpen ? "Hide" : "Settings"}
      </button>

      <div
        id="sim-menu-panel"
        className="settings-menu"
        role="dialog"
        aria-label="Settings"
      >
        <div className="settings-title">
          Settings (Left Click On Simulation To Affect)
        </div>
        <div className="settings-section">
          <div className="settings-row">
            <label htmlFor="density">Rest Density</label>
            <div className="settings-control">
              <input
                id="density"
                type="range"
                min="0"
                max="30"
                step="0.01"
                value={restDensity}
                onChange={(event) => setRestDensity(Number(event.target.value))}
              />
              <span className="settings-value">{restDensity.toFixed(2)}</span>
            </div>
          </div>

          <div className="settings-row">
            <label htmlFor="gas-const">Gas Constant</label>
            <div className="settings-control">
              <input
                id="gas-const"
                type="range"
                min="1"
                max="180"
                step="0.01"
                value={gasConst}
                onChange={(event) => setGasConst(Number(event.target.value))}
              />
              <span className="settings-value">{gasConst.toFixed(2)}</span>
            </div>
          </div>

          <div className="settings-row">
            <label htmlFor="gravity">Gravity</label>
            <div className="settings-control">
              <input
                id="gravity"
                type="range"
                min="0"
                max="50"
                step="0.01"
                value={gravity}
                onChange={(event) => setGravity(Number(event.target.value))}
              />
              <span className="settings-value">{gravity.toFixed(2)}</span>
            </div>
          </div>

          <div className="settings-row">
            <label htmlFor="viscocity">Viscosity</label>
            <div className="settings-control">
              <input
                id="viscocity"
                type="range"
                min="14"
                max="100"
                step="0.01"
                value={viscosity}
                onChange={(event) => setViscosity(Number(event.target.value))}
              />
              <span className="settings-value">{viscosity.toFixed(2)}</span>
            </div>
          </div>

          <div className="settings-row">
            <label htmlFor="particles">Particle Count</label>
            <div className="settings-control">
              <input
                id="particles"
                type="range"
                min="0"
                max="10000"
                step="1"
                value={particleCount}
                onChange={(event) => setParticleCount(Number(event.target.value))}
              />
              <span className="settings-value">
                {particleCount.toLocaleString()}
              </span>
            </div>
          </div>

          <div className="settings-row settings-row--toggle">
            <label htmlFor="mouse-behaviour">Toggle Mouse Behaviour</label>
            <label className="toggle-switch">
              <input
                id="mouse-behaviour"
                type="checkbox"
                checked={mouseEffect}
                onChange={(event) => setMouseEffect(event.target.checked)}
              />
              <span className="toggle-slider" />
            </label>
          </div>
        </div>
      </div>
    </div>
  );
}

export default SimMenu;