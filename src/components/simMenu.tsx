import "../css/simMenu.css";

function SimMenu() {
  return (
    <div className="settings-menu" role="dialog" aria-label="Settings">
      <div className="settings-title">
        Settings (Left Click On Simulation To Affect)
      </div>
      <div className="settings-section">
        <div className="settings-row">
          <label htmlFor="density">Rest Density</label>
          <input id="density" type="range" min="0" max="30" step="0.01" defaultValue="4" />
        </div>

        <div className="settings-row">
          <label htmlFor="gas-const">Gas Constant</label>
          <input id="gas-const" type="range" min="1" max="180" step="0.01" defaultValue="80" />
        </div>

        <div className="settings-row">
          <label htmlFor="gravity">Gravity</label>
          <input id="gravity" type="range" min="0" max="50" step="0.01" defaultValue="9" />
        </div>

        <div className="settings-row">
          <label htmlFor="viscocity">Viscosity</label>
          <input id="viscocity" type="range" min="14" max="100" step="0.01" defaultValue="80" />
        </div>

        <div className="settings-row">
          <label htmlFor="particles">Particle Count</label>
          <input id="particles" type="range" min="0" max="10000" step="1" defaultValue="3000" />
        </div>
      </div>
    </div>
  );
}

export default SimMenu;