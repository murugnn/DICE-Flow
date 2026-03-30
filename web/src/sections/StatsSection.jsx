import { useEffect, useRef, useState } from 'react';
import { motion, useInView } from 'framer-motion';
import { dashboardStats } from '../data/mockData';

// Animated counter hook
const useCounter = (end, duration = 2000, startOnView = false) => {
  const [count, setCount] = useState(0);
  const ref = useRef(null);
  const isInView = useInView(ref, { once: true, margin: "-100px" });
  const hasStarted = useRef(false);

  useEffect(() => {
    if ((startOnView && !isInView) || hasStarted.current) return;
    
    hasStarted.current = true;
    let startTime = null;
    const startValue = 0;

    const animate = (timestamp) => {
      if (!startTime) startTime = timestamp;
      const progress = Math.min((timestamp - startTime) / duration, 1);
      
      // Easing function (ease-out-cubic)
      const easeOut = 1 - Math.pow(1 - progress, 3);
      setCount(Math.floor(startValue + (end - startValue) * easeOut));

      if (progress < 1) {
        requestAnimationFrame(animate);
      }
    };

    requestAnimationFrame(animate);
  }, [end, duration, isInView, startOnView]);

  return { count, ref };
};

const StatCard = ({ value, label, suffix = '', prefix = '', delay = 0 }) => {
  const { count, ref } = useCounter(value, 2000, true);

  return (
    <motion.div
      ref={ref}
      className="glass-card p-6 flex flex-col items-center justify-center text-center"
      initial={{ opacity: 0, y: 40, scale: 0.98 }}
      whileInView={{ opacity: 1, y: 0, scale: 1 }}
      viewport={{ once: true, margin: "-50px" }}
      transition={{ duration: 0.6, delay, ease: [0.25, 0.46, 0.45, 0.94] }}
    >
      <div className="font-mono-data text-[clamp(32px,4vw,48px)] font-bold text-[#00F0FF] text-glow-cyan">
        {prefix}{count.toLocaleString()}{suffix}
      </div>
      <div className="mt-2 text-sm text-[#A6AEB8] uppercase tracking-wider">
        {label}
      </div>
    </motion.div>
  );
};

const StatsSection = () => {
  const sectionRef = useRef(null);

  const stats = [
    { value: dashboardStats.edgesProcessed, label: 'Network edges analyzed', suffix: '+' },
    { value: dashboardStats.nodesAnalyzed, label: 'Genes scored' },
    { value: 10, label: 'Rewiring candidates', prefix: 'Top ' }
  ];

  const features = [
    'Differential Betweenness & Eigenvector',
    'Ensemble ranking with outlier rejection',
    'Export-ready gene lists & subnetworks'
  ];

  return (
    <section 
      ref={sectionRef}
      className="relative w-full min-h-screen bg-[#0A0A0A] py-20 px-4 sm:px-8 lg:px-16"
    >
      {/* Subtle mesh gradient */}
      <div className="absolute inset-0 gradient-mesh pointer-events-none opacity-50" />

      <div className="relative z-10 max-w-7xl mx-auto">
        {/* Stats Row */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-20">
          {stats.map((stat, i) => (
            <StatCard
              key={stat.label}
              value={stat.value}
              label={stat.label}
              suffix={stat.suffix}
              prefix={stat.prefix}
              delay={i * 0.1}
            />
          ))}
        </div>

        {/* Content Grid */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 lg:gap-20">
          {/* Left Column - Title & Description */}
          <motion.div
            initial={{ opacity: 0, x: -40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true, margin: "-100px" }}
            transition={{ duration: 0.7, ease: [0.25, 0.46, 0.45, 0.94] }}
          >
            <h2 className="font-display text-[clamp(28px,3.5vw,44px)] font-bold text-[#F2F4F8] leading-tight">
              What DiCE-Duan measures
            </h2>
            <p className="mt-6 text-[#A6AEB8] leading-relaxed text-base lg:text-lg">
              We compare Normal and Tumor networks to find genes that change their influence—
              rewiring events that can drive progression. DiCE-Duan combines differential 
              centrality with an ensemble of topological signatures to rank the most 
              actionable switches.
            </p>
          </motion.div>

          {/* Right Column - Features */}
          <motion.div
            initial={{ opacity: 0, x: 40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true, margin: "-100px" }}
            transition={{ duration: 0.7, delay: 0.1, ease: [0.25, 0.46, 0.45, 0.94] }}
            className="space-y-4"
          >
            {features.map((feature, i) => (
              <motion.div
                key={feature}
                className="flex items-center gap-4 p-4 rounded-xl border border-white/[0.06] bg-white/[0.02]"
                initial={{ opacity: 0, y: 20 }}
                whileInView={{ opacity: 1, y: 0 }}
                viewport={{ once: true }}
                transition={{ duration: 0.5, delay: 0.2 + i * 0.1 }}
              >
                <div className="w-8 h-8 rounded-lg bg-[#00F0FF]/10 flex items-center justify-center flex-shrink-0">
                  <svg className="w-4 h-4 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                  </svg>
                </div>
                <span className="text-[#F2F4F8] text-sm lg:text-base">{feature}</span>
              </motion.div>
            ))}

            {/* Pipeline Diagram Card */}
            <motion.div
              className="mt-8 glass-card p-6"
              initial={{ opacity: 0, y: 30 }}
              whileInView={{ opacity: 1, y: 0 }}
              viewport={{ once: true }}
              transition={{ duration: 0.6, delay: 0.5 }}
            >
              <div className="flex items-center justify-between gap-2">
                {['Input', 'Network', 'Score', 'Rank'].map((step, i) => (
                  <div key={step} className="flex items-center gap-2">
                    <div className="flex flex-col items-center">
                      <div 
                        className="w-10 h-10 rounded-full flex items-center justify-center text-xs font-mono-data font-bold"
                        style={{
                          background: i === 0 ? 'rgba(0, 240, 255, 0.2)' : 'rgba(255, 255, 255, 0.05)',
                          border: `1px solid ${i === 0 ? 'rgba(0, 240, 255, 0.4)' : 'rgba(255, 255, 255, 0.1)'}`,
                          color: i === 0 ? '#00F0FF' : '#A6AEB8'
                        }}
                      >
                        {i + 1}
                      </div>
                      <span className="mt-2 text-[10px] text-[#A6AEB8] uppercase tracking-wider">{step}</span>
                    </div>
                    {i < 3 && (
                      <svg className="w-4 h-4 text-[#A6AEB8]/50" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 5l7 7-7 7" />
                      </svg>
                    )}
                  </div>
                ))}
              </div>
            </motion.div>
          </motion.div>
        </div>
      </div>
    </section>
  );
};

export default StatsSection;
