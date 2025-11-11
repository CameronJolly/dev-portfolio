import { useState } from "react";
import { BrowserRouter, Routes, Route } from "react-router-dom";
import "./css/app.css";
import MainLayout from "./pages/mainLayout";
import Dashboard from "./pages/homePage";
import Contact from "./pages/contact";
import Portfolio from "./pages/portfolio";
import IntroVideo from "./components/intro";
import Resume from "./pages/resume";
import ParticleSim from "./pages/particleSim";

function App() {
  const [showIntro, setShowIntro] = useState(true);

  return (
    <>
      <BrowserRouter>
        <Routes>
          <Route element={<MainLayout />}>
            <Route path="/" element={<Dashboard />} />
            <Route path="/Portfolio" element={<Portfolio />} />
            <Route path="/portfolio" element={<Portfolio/>} />
            <Route path="/Contact" element={<Contact />} />
            <Route path="/contact" element={<Contact/>} />
            <Route path="/Resume" element={<Resume />}/>
            <Route path="/resume" element={<Resume/>}/>
            <Route path="/ParticleSim" element={<ParticleSim/>}/>
            <Route path="/Particlesim" element={<ParticleSim/>}/>
            <Route path="/particleSim" element={<ParticleSim/>}/>
            <Route path="/particlesim" element={<ParticleSim/>}/>
          </Route>
        </Routes>
      </BrowserRouter>
      {showIntro && <IntroVideo onFinish={() => setShowIntro(false)} />}
    </>
  );
}

export default App;
