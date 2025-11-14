import { useEffect, useRef } from "react";
import { useLocation } from "react-router-dom";
import "../css/background.css";

type Particle = {
  x: number;
  y: number;
  vx: number;
  vy: number;
  size: number;
};

export default function BackgroundCanvas() {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const location = useLocation();

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const dpr = Math.min(window.devicePixelRatio || 1, 2);

    const measureAndResizeCanvas = () => {
      const prevDisplay = canvas.style.display;
      canvas.style.display = "none";

      const docEl = document.documentElement;
      const body = document.body;

      const contentWidth = Math.max(
        docEl.scrollWidth,
        docEl.clientWidth,
        body?.scrollWidth || 0
      );

      const contentHeight = Math.max(
        docEl.scrollHeight,
        docEl.clientHeight,
        body?.scrollHeight || 0,
        window.innerHeight
      );

      canvas.style.display = prevDisplay || "";

      canvas.style.width = `${contentWidth}px`;
      canvas.style.height = `${contentHeight}px`;

      canvas.width = Math.floor(contentWidth * dpr);
      canvas.height = Math.floor(contentHeight * dpr);

      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    };

    measureAndResizeCanvas();

    const onWindowResize = () => measureAndResizeCanvas();
    window.addEventListener("resize", onWindowResize);

    const ro = new ResizeObserver(() => measureAndResizeCanvas());
    ro.observe(document.body);

    const rafMeasure = requestAnimationFrame(measureAndResizeCanvas);

    const mobileCheck = function() {
      let check = false;
      (function(a){if(/(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino/i.test(a)||/1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas\-|your|zeto|zte\-/i.test(a.substr(0,4))) check = true;})(navigator.userAgent||navigator.vendor);
      return check;
    };

    const mobileStatus = mobileCheck();

    let particleAmount = 200;

    if (mobileStatus === true) {
      particleAmount = 50;
    }

    const particles: Particle[] = Array.from({ length: particleAmount }, () => ({
      x: Math.random() * (canvas.width / dpr),
      y: Math.random() * (canvas.height / dpr),
      vx: (Math.random() - 0.5) * 0.2,
      vy: (Math.random() - 0.5) * 0.2,
      size: Math.random() * 1.7 + 1,
    }));

    const mouse = { x: 0, y: 0 };

    const handleMouseMove = (e: MouseEvent) => {
      const rect = canvas.getBoundingClientRect();
      mouse.x = e.clientX - rect.left;
      mouse.y = e.clientY - rect.top;
    };
    window.addEventListener("mousemove", handleMouseMove);

    let rafId = 0;
    let running = true;

    const draw = () => {
      if (!running) return;

      ctx.fillStyle = "white";
      ctx.fillRect(0, 0, canvas.width / dpr, canvas.height / dpr);

      for (let p of particles) {
        const dx = p.x - mouse.x;
        const dy = p.y - mouse.y;

        if (Math.abs(dx) < 200 && Math.abs(dy) < 200) {
          p.x += dx * 0.0009;
          p.y += dy * 0.0009;
        }

        p.x += p.vx;
        p.y += p.vy;

        if (p.x < 0 || p.x > canvas.width / dpr) p.vx *= -1;
        if (p.y < 0 || p.y > canvas.height / dpr) p.vy *= -1;

        for (let i of particles) {
          if (Math.abs(p.x - i.x) < 100 && Math.abs(p.y - i.y) < 100) {
            ctx.beginPath();
            ctx.moveTo(p.x, p.y);
            ctx.lineTo(i.x, i.y);
            ctx.lineWidth = 0.8;
            ctx.strokeStyle = "rgba(26, 26, 26, 0.53)";
            ctx.stroke();
          }
        }

        ctx.beginPath();
        ctx.arc(p.x, p.y, p.size, 0, Math.PI * 2);
        ctx.fillStyle = "rgba(40, 40, 40, 1)";
        ctx.fill();
      }

      rafId = requestAnimationFrame(draw);
    };

    draw();

    
    return () => {
      running = false;
      cancelAnimationFrame(rafId);
      cancelAnimationFrame(rafMeasure);
      window.removeEventListener("resize", onWindowResize);
      window.removeEventListener("mousemove", handleMouseMove);
      ro.disconnect();
    };
  }, [location]); 

  return <canvas ref={canvasRef} className="background-canvas" />;
}
