import { NavLink } from "react-router-dom";
import "../css/homepage.css";
import { motion } from "framer-motion";
import profilePic from "../assets/Profile-Picture.png";
import { useEffect, useState } from "react";

// Portfolio images for preloading
import BudgetingPic from "../assets/BudgetingTool.jpg";
import CameronPic from "../assets/cameronjolly.com.jpg";
import FluidSim from "../assets/FluidSim.jpg";
import WebSim from "../assets/webparticlesim.png";

function Homepage() {
  const [hasVisitedBefore, setHasVisitedBefore] = useState(
    sessionStorage.getItem("hasVisited") === "true",
  );

  // Mark user visit and preload images
  useEffect(() => {
    // Mark that user has visited
    if (!hasVisitedBefore) {
      sessionStorage.setItem("hasVisited", "true");
      setHasVisitedBefore(true);
    }
  }, [hasVisitedBefore]);

  // Preload portfolio images silently
  useEffect(() => {
    const imagesToPreload = [BudgetingPic, CameronPic, FluidSim, WebSim];
    imagesToPreload.forEach((src) => {
      const img = new Image();
      img.src = src;
    });
  }, []);

  const containerVariants = {
    hidden: { opacity: 0 },
    visible: {
      opacity: 1,
      transition: {
        staggerChildren: 0.2,
        delayChildren: hasVisitedBefore ? 0 : 1.3,
      },
    },
  };

  const itemVariants = {
    hidden: { y: 20, opacity: 0 },
    visible: {
      y: 0,
      opacity: 1,
      transition: {
        duration: 0.8,
        ease: [0.25, 0.1, 0.25, 1] as const,
      },
    },
  };

  return (
    <div className="home-content">
      <div className="backdrop"></div>
      <motion.div
        className="hero-container"
        variants={containerVariants}
        initial="hidden"
        animate="visible"
      >
        {/* Hero Section */}
        <motion.div className="hero-section" variants={itemVariants}>
          <div className="hero-text">
            <motion.div
              className="greeting"
              initial={{ x: -50, opacity: 0 }}
              animate={{ x: 0, opacity: 1 }}
              transition={{ duration: 0.6, delay: hasVisitedBefore ? 0 : 1.5 }}
            >
              Hi, I'm
            </motion.div>
            <motion.h1
              className="hero-name"
              initial={{ x: -50, opacity: 0 }}
              animate={{ x: 0, opacity: 1 }}
              transition={{
                duration: 0.8,
                delay: hasVisitedBefore ? 0.1 : 1.7,
              }}
            >
              Cameron Jolly
            </motion.h1>
            <motion.p
              className="hero-title"
              initial={{ x: -50, opacity: 0 }}
              animate={{ x: 0, opacity: 1 }}
              transition={{
                duration: 0.8,
                delay: hasVisitedBefore ? 0.2 : 1.9,
              }}
            >
              Full Stack Developer & Creative Technologist
            </motion.p>
            <motion.div
              className="hero-cta"
              initial={{ y: 20, opacity: 0 }}
              animate={{ y: 0, opacity: 1 }}
              transition={{
                duration: 0.8,
                delay: hasVisitedBefore ? 0.3 : 2.1,
              }}
            >
              <NavLink to="/portfolio" className="cta-button primary">
                View My Work
              </NavLink>
              <NavLink to="/contact" className="cta-button secondary">
                Get in Touch
              </NavLink>
            </motion.div>
          </div>

          <motion.div
            className="hero-image"
            initial={{ scale: 0.8, opacity: 0 }}
            animate={{ scale: 1, opacity: 1 }}
            transition={{ duration: 1, delay: hasVisitedBefore ? 0.2 : 1.8 }}
          >
            <div className="image-wrapper">
              <img src={profilePic} alt="Cameron Jolly" />
            </div>
          </motion.div>
        </motion.div>

        {/* About Section */}
        <motion.div className="about-section" variants={itemVariants}>
          <div className="section-header">
            <h2>About Me</h2>
            <div className="header-line"></div>
          </div>
          <div className="about-content">
            <p className="about-text">
              I'm a Full Stack Web Development student and creative technologist
              passionate about building responsive, interactive web experiences.
              I specialize in modern technologies like <strong>React</strong>,{" "}
              <strong>TypeScript</strong>,<strong>Python</strong>, and{" "}
              <strong>C++</strong>.
            </p>
            <p className="about-text">
              With a strong foundation in problem-solving and a collaborative
              mindset, I combine technical expertise with a background in design
              to create engaging digital experiences that make an impact.
            </p>
          </div>
        </motion.div>

        {/* Skills Section */}
        <motion.div className="skills-section" variants={itemVariants}>
          <div className="section-header">
            <h2>Technical Skills</h2>
            <div className="header-line"></div>
          </div>
          <div className="skills-grid">
            <div className="skill-card">
              <h3>Frontend</h3>
              <p>React, TypeScript, JavaScript, HTML, CSS</p>
            </div>
            <div className="skill-card">
              <h3>Backend</h3>
              <p>Python, C++, Node.js, APIs</p>
            </div>
            <div className="skill-card">
              <h3>Design</h3>
              <p>UI/UX, Responsive Design, Animations</p>
            </div>
            <div className="skill-card">
              <h3>Tools</h3>
              <p>Git, VS Code, Vite, Framer Motion</p>
            </div>
          </div>
        </motion.div>

        {/* Contact CTA Section */}
        <motion.div className="contact-cta-section" variants={itemVariants}>
          <h2>Let's Build Something Together</h2>
          <p>I'm always open to discussing new projects and opportunities.</p>
          <div className="contact-links">
            <a href="mailto:cameron.l.jolly@gmail.com" className="contact-link">
              cameron.l.jolly@gmail.com
            </a>
            <NavLink to="/resume" className="contact-link">
              View Resume
            </NavLink>
          </div>
        </motion.div>
      </motion.div>
    </div>
  );
}

export default Homepage;
