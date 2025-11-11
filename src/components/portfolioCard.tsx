import "../css/portfolioCard.css"
import { NavLink } from 'react-router-dom';


function PortfolioCard({description='',imageLink='', title='', projectLink = ''}) {
    const isExternal = projectLink.startsWith('http');

    return(
        <div className="portfolio-card">
            <div className="card-image">
                <img src={imageLink} alt={title} />
            </div>
            <div className="card-info">
                <h1 className="card-title">{title}</h1>
                <p className="card-description">{description}</p>
            </div>
            <div className="card-action">
                {isExternal ? (
                    <a target="_blank" rel="noopener noreferrer" href={projectLink}><button>→</button></a>
                ) : (
                    <NavLink to={projectLink}><button>→</button></NavLink>
                )}
                <span className="check-out-text">Check it Out</span>
            </div>
        </div>
    )    
}

export default PortfolioCard