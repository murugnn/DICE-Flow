import { useEffect, useRef, useState } from 'react';
import { motion } from 'framer-motion';

const HeroSection = ({ onExploreClick }) => {
  const canvasRef = useRef(null);
  const animationRef = useRef(null);
  const [isLoaded, setIsLoaded] = useState(false);

  // Network background animation
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    let width = canvas.width = canvas.offsetWidth;
    let height = canvas.height = canvas.offsetHeight;

    // Network nodes
    const nodes = [];
    const nodeCount = 50;
    const connectionDistance = 120;

    for (let i = 0; i < nodeCount; i++) {
      nodes.push({
        x: Math.random() * width,
        y: Math.random() * height,
        vx: (Math.random() - 0.5) * 0.3,
        vy: (Math.random() - 0.5) * 0.3,
        size: Math.random() * 3 + 2,
        pulsePhase: Math.random() * Math.PI * 2
      });
    }

    let frameCount = 0;
    const animate = () => {
      frameCount++;
      // Render every 2nd frame for performance
      if (frameCount % 2 === 0) {
        ctx.fillStyle = 'rgba(10, 10, 10, 0.15)';
        ctx.fillRect(0, 0, width, height);

        // Update and draw nodes
        nodes.forEach((node, i) => {
          // Update position
          node.x += node.vx;
          node.y += node.vy;

          // Bounce off edges
          if (node.x < 0 || node.x > width) node.vx *= -1;
          if (node.y < 0 || node.y > height) node.vy *= -1;

          // Keep in bounds
          node.x = Math.max(5, Math.min(width - 5, node.x));
          node.y = Math.max(5, Math.min(height - 5, node.y));

          // Draw connections (only check every 5th node for performance)
          if (i % 5 === 0) {
            nodes.forEach((other, j) => {
              if (i >= j) return;
              const dx = node.x - other.x;
              const dy = node.y - other.y;
              const dist = Math.sqrt(dx * dx + dy * dy);

              if (dist < connectionDistance) {
                const opacity = (1 - dist / connectionDistance) * 0.3;
                ctx.beginPath();
                ctx.moveTo(node.x, node.y);
                ctx.lineTo(other.x, other.y);
                ctx.strokeStyle = `rgba(0, 240, 255, ${opacity})`;
                ctx.lineWidth = 0.5;
                ctx.stroke();
              }
            });
          }

          // Draw node with pulse
          const pulse = Math.sin(Date.now() * 0.002 + node.pulsePhase) * 0.3 + 0.7;
          ctx.beginPath();
          ctx.arc(node.x, node.y, node.size * pulse, 0, Math.PI * 2);
          ctx.fillStyle = `rgba(0, 240, 255, ${0.4 + pulse * 0.4})`;
          ctx.fill();

          // Glow effect for some nodes
          if (i % 7 === 0) {
            ctx.beginPath();
            ctx.arc(node.x, node.y, node.size * 3, 0, Math.PI * 2);
            const gradient = ctx.createRadialGradient(
              node.x, node.y, 0,
              node.x, node.y, node.size * 3
            );
            gradient.addColorStop(0, 'rgba(0, 240, 255, 0.2)');
            gradient.addColorStop(1, 'rgba(0, 240, 255, 0)');
            ctx.fillStyle = gradient;
            ctx.fill();
          }
        });
      }

      animationRef.current = requestAnimationFrame(animate);
    };

    animate();
    setIsLoaded(true);

    const handleResize = () => {
      width = canvas.width = canvas.offsetWidth;
      height = canvas.height = canvas.offsetHeight;
    };

    window.addEventListener('resize', handleResize);

    return () => {
      cancelAnimationFrame(animationRef.current);
      window.removeEventListener('resize', handleResize);
    };
  }, []);

  // Text animation variants
  const containerVariants = {
    hidden: { opacity: 0 },
    visible: {
      opacity: 1,
      transition: {
        staggerChildren: 0.08,
        delayChildren: 0.3
      }
    }
  };

  const charVariants = {
    hidden: { 
      opacity: 0, 
      y: 24,
      filter: 'blur(8px)'
    },
    visible: { 
      opacity: 1, 
      y: 0,
      filter: 'blur(0px)',
      transition: {
        duration: 0.5,
        ease: [0.25, 0.46, 0.45, 0.94]
      }
    }
  };

  const title = "DiCE-Duan Analyzer";

  return (
    <section className="relative w-full h-screen overflow-hidden bg-[#0A0A0A]">
      {/* Background Network Canvas */}
      <canvas
        ref={canvasRef}
        className="absolute inset-0 w-full h-full"
        style={{ opacity: isLoaded ? 1 : 0, transition: 'opacity 0.8s ease' }}
      />

      {/* Cyan Bloom Effect */}
      <div 
        className="absolute left-1/2 top-[38%] -translate-x-1/2 -translate-y-1/2 
                   w-[90vw] h-[70vh] pointer-events-none"
        style={{
          background: 'radial-gradient(ellipse at center, rgba(0, 240, 255, 0.12) 0%, transparent 60%)',
          filter: 'blur(60px)'
        }}
      />

      {/* Content */}
      <div className="relative z-10 flex flex-col items-center justify-center h-full px-4">
        {/* Main Title */}
        <motion.h1 
          className="font-display text-[clamp(44px,8vw,72px)] font-bold text-[#F2F4F8] text-center tracking-tight"
          variants={containerVariants}
          initial="hidden"
          animate="visible"
        >
          {title.split('').map((char, i) => (
            <motion.span
              key={i}
              variants={charVariants}
              className={char === ' ' ? 'inline' : 'inline-block'}
              style={{ 
                textShadow: '0 0 40px rgba(0, 240, 255, 0.3)',
                marginRight: char === ' ' ? '0.3em' : '0'
              }}
            >
              {char}
            </motion.span>
          ))}
        </motion.h1>

        {/* Subtitle */}
        <motion.p
          className="mt-6 text-[clamp(14px,2vw,18px)] text-[#A6AEB8] text-center max-w-2xl leading-relaxed"
          initial={{ opacity: 0, y: 14 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.8 }}
        >
          Differential Centrality & Ensemble rewiring analysis for cancer networks
        </motion.p>

        {/* CTA Button */}
        <motion.button
          onClick={onExploreClick}
          className="mt-10 btn-glass flex items-center gap-3 group"
          initial={{ opacity: 0, scale: 0.96 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.4, delay: 1 }}
          whileHover={{ scale: 1.02 }}
          whileTap={{ scale: 0.98 }}
        >
          <span>Explore the network</span>
          <svg 
            className="w-4 h-4 transition-transform duration-300 group-hover:translate-x-1" 
            fill="none" 
            viewBox="0 0 24 24" 
            stroke="currentColor"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 8l4 4m0 0l-4 4m4-4H3" />
          </svg>
        </motion.button>

        {/* Scroll Hint */}
        <motion.div
          className="absolute bottom-8 left-1/2 -translate-x-1/2 flex flex-col items-center gap-2"
          initial={{ opacity: 0 }}
          animate={{ opacity: 0.6 }}
          transition={{ duration: 0.5, delay: 1.5 }}
        >
          <span className="text-xs text-[#A6AEB8] uppercase tracking-widest">Scroll to analyze</span>
          <motion.div
            className="w-5 h-8 rounded-full border border-[#A6AEB8]/30 flex justify-center pt-2"
            animate={{ y: [0, 4, 0] }}
            transition={{ duration: 1.5, repeat: Infinity, ease: "easeInOut" }}
          >
            <div className="w-1 h-2 bg-[#00F0FF] rounded-full" />
          </motion.div>
        </motion.div>
      </div>

      {/* Top Navigation */}
      <motion.nav
        className="absolute top-0 left-0 right-0 z-20 flex items-center justify-between px-8 py-6"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5, delay: 0.2 }}
      >
        <div className="font-display text-xl font-bold text-[#F2F4F8]">
          DiCE-Duan
        </div>
        <div className="hidden md:flex items-center gap-8">
          {['Explorer', 'Leaderboard', 'Method', 'Contact'].map((item, i) => (
            <a
              key={item}
              href={`#${item.toLowerCase()}`}
              className="text-sm text-[#A6AEB8] hover:text-[#00F0FF] transition-colors duration-200"
            >
              {item}
            </a>
          ))}
        </div>
      </motion.nav>
    </section>
  );
};

export default HeroSection;
