import "../css/portfolio.css"
import PortfolioCard from "../components/portfolioCard"
import BudgetingPic from '../assets/BudgetingTool.jpg'
import CameronPic from '../assets/cameronjolly.com.jpg'
import FluidSim from '../assets/FluidSim.jpg'
import WebSim from '../assets/webparticlesim.png'

function Portfolio() {
    return (
        <>
            <div className="backdrop">
            </div>

            <h1 className="port-header">My Work</h1>

            <div className="port-card-holder">
                <PortfolioCard projectLink="/ParticleSim" title="Interactive Web SPH Fluid Sim" imageLink={WebSim} description="An interactive SPH particle simulation created using ThreeJS for rendering and a custom physics engine     (For the best experience use Chrome)"/>
                <PortfolioCard projectLink="https://cameronjolly.com/" title="Portfolio Website" imageLink={CameronPic} description="A portfolio website for displaying photography, this includes a product page and contact form made in vanilla HTML/CSS"/>
                <PortfolioCard projectLink="https://github.com/Majorkamo02/Python-Budgeting-Program" title="Budgeting Tool" imageLink={BudgetingPic} description="A local budgeting program utilizing encrypted logins, a SQLite database, and dynamic data visualizations"/>
                <PortfolioCard projectLink="https://github.com/Majorkamo02/Cpp-CUDA-Simulation" title="SPH Fluid Simulation" imageLink={FluidSim} description="A parallelized SPH Fluid simulation in C++ using the CUDA architecture for performance and SFML for the particle visualization"/>
            </div>
        </>
    )
}

export default Portfolio