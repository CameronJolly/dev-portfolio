import { useState, useEffect } from "react";
import { BrowserRouter, Routes, Route, useLocation } from "react-router-dom";
import "./css/app.css";
import MainLayout from "./pages/mainLayout";
import Dashboard from "./pages/homePage";
import Contact from "./pages/contact";
import Portfolio from "./pages/portfolio";
import IntroVideo from "./components/intro";
import Resume from "./pages/resume";
import ParticleSim from "./pages/particleSim";

function AppContent() {
  const location = useLocation();
  const [showIntro, setShowIntro] = useState(false);

  useEffect(() => {
    // Only show intro if on homepage and hasn't been shown this session
    const hasSeenIntro = sessionStorage.getItem("hasSeenIntro");
    if (location.pathname === "/" && !hasSeenIntro) {
      setShowIntro(true);
      sessionStorage.setItem("hasSeenIntro", "true");
    }
  }, [location.pathname]);

  return (
    <>
      <Routes>
        <Route element={<MainLayout />}>
          <Route path="/" element={<Dashboard />} />
          <Route path="/Portfolio" element={<Portfolio />} />
          <Route path="/portfolio" element={<Portfolio />} />
          <Route path="/Contact" element={<Contact />} />
          <Route path="/contact" element={<Contact />} />
          <Route path="/Resume" element={<Resume />} />
          <Route path="/resume" element={<Resume />} />
          <Route path="/ParticleSim" element={<ParticleSim />} />
          <Route path="/Particlesim" element={<ParticleSim />} />
          <Route path="/particleSim" element={<ParticleSim />} />
          <Route path="/particlesim" element={<ParticleSim />} />
        </Route>
      </Routes>
      {showIntro && <IntroVideo onFinish={() => setShowIntro(false)} />}
    </>
  );
}

function App() {
  return (
    <BrowserRouter>
      <AppContent />
    </BrowserRouter>
  );
}

export default App;
