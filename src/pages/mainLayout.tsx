import { Outlet, useLocation } from "react-router-dom";
import Header from "../components/header";
import Footer from "../components/footer";
import BackgroundCanvas from "../components/background";

export default function MainLayout() {
  const { pathname } = useLocation();
  const hideFooter = pathname.toLowerCase().includes("particlesim");

  return (
    <>     
      <BackgroundCanvas /> 
      {!hideFooter && <Header />}
      <main>
        <Outlet /> {}
      </main>
      {!hideFooter && <Footer />}
    </>
  );
}