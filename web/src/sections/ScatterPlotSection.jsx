import { useState } from 'react';
import { motion } from 'framer-motion';
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ZAxis,
  ReferenceLine
} from 'recharts';
import { scatterData } from '../data/mockData';

const CustomTooltip = ({ active, payload }) => {
  if (active && payload && payload.length) {
    const data = payload[0].payload;
    return (
      <div className="glass-card-sm p-4 border border-[#00F0FF]/30">
        <p className="font-mono-data font-bold text-[#F2F4F8] mb-2">{data.gene}</p>
        <div className="space-y-1 text-xs">
          <p className="text-[#A6AEB8]">
            Delta Betweenness: <span className="text-[#00F0FF]">{data.deltaBetweenness.toFixed(3)}</span>
          </p>
          <p className="text-[#A6AEB8]">
            Delta Eigenvector: <span className="text-[#00F0FF]">{data.deltaEigenvector.toFixed(3)}</span>
          </p>
          <p className="text-[#A6AEB8]">
            Ensemble Score: <span className="text-[#FFD166]">{data.ensembleScore.toFixed(2)}</span>
          </p>
        </div>
      </div>
    );
  }
  return null;
};

const ScatterPlotSection = () => {
  const [hoveredPoint, setHoveredPoint] = useState(null);

  // Separate data by category
  const topGenes = scatterData.filter(d => d.category === 'top');
  const candidateGenes = scatterData.filter(d => d.category === 'candidate');
  const otherGenes = scatterData.filter(d => d.category === 'other');

  return (
    <section className="relative w-full min-h-screen bg-[#0A0A0A] py-20 px-4 sm:px-8 lg:px-16">
      {/* Background gradient */}
      <div className="absolute inset-0 gradient-mesh pointer-events-none opacity-30" />

      <div className="relative z-10 max-w-7xl mx-auto">
        <div className="grid grid-cols-1 lg:grid-cols-12 gap-8 lg:gap-12">
          {/* Left Column - Description */}
          <motion.div
            className="lg:col-span-4"
            initial={{ opacity: 0, x: -40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true, margin: "-100px" }}
            transition={{ duration: 0.6 }}
          >
            <h2 className="font-display text-[clamp(24px,3vw,36px)] font-bold text-[#F2F4F8] leading-tight">
              Centrality shift
              <span className="block text-[#A6AEB8]">(Normal vs Tumor)</span>
            </h2>
            
            <p className="mt-6 text-[#A6AEB8] leading-relaxed">
              Each point represents a gene. Axes show the change in Betweenness and 
              Eigenvector centrality between Normal and Tumor conditions. Genes in the 
              upper-right quadrant show increased influence in tumor networks.
            </p>

            {/* Legend */}
            <div className="mt-8 space-y-3">
              <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider">
                Legend
              </h4>
              
              <div className="flex items-center gap-3">
                <div className="w-4 h-4 rounded-full bg-[#FFD166] shadow-lg"
                     style={{ boxShadow: '0 0 12px rgba(255, 209, 102, 0.5)' }} />
                <span className="text-sm text-[#A6AEB8]">Top rewired genes</span>
              </div>
              
              <div className="flex items-center gap-3">
                <div className="w-3 h-3 rounded-full bg-[#00F0FF]"
                     style={{ boxShadow: '0 0 8px rgba(0, 240, 255, 0.4)' }} />
                <span className="text-sm text-[#A6AEB8]">Candidate genes</span>
              </div>
              
              <div className="flex items-center gap-3">
                <div className="w-2 h-2 rounded-full bg-[#A6AEB8]/50" />
                <span className="text-sm text-[#A6AEB8]">Other genes</span>
              </div>
            </div>

            {/* Stats */}
            <div className="mt-8 grid grid-cols-2 gap-4">
              <div className="glass-card-sm p-4">
                <div className="text-xs text-[#A6AEB8] uppercase tracking-wider mb-1">
                  Positive Shift
                </div>
                <div className="font-mono-data text-xl font-bold text-[#00F0FF]">
                  {scatterData.filter(d => d.deltaBetweenness > 0 && d.deltaEigenvector > 0).length}
                </div>
                <div className="text-xs text-[#A6AEB8]/70 mt-1">genes</div>
              </div>
              
              <div className="glass-card-sm p-4">
                <div className="text-xs text-[#A6AEB8] uppercase tracking-wider mb-1">
                  Max Delta
                </div>
                <div className="font-mono-data text-xl font-bold text-[#FFD166]">
                  {Math.max(...scatterData.map(d => d.deltaBetweenness)).toFixed(2)}
                </div>
                <div className="text-xs text-[#A6AEB8]/70 mt-1">betweenness</div>
              </div>
            </div>
          </motion.div>

          {/* Right Column - Scatter Plot */}
          <motion.div
            className="lg:col-span-8"
            initial={{ opacity: 0, x: 40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true, margin: "-100px" }}
            transition={{ duration: 0.6, delay: 0.1 }}
          >
            <div className="glass-card p-6" style={{ height: '500px' }}>
              <ResponsiveContainer width="100%" height="100%">
                <ScatterChart
                  margin={{ top: 20, right: 30, bottom: 40, left: 50 }}
                >
                  <CartesianGrid 
                    strokeDasharray="3 3" 
                    stroke="rgba(255, 255, 255, 0.06)"
                  />
                  
                  <XAxis
                    type="number"
                    dataKey="deltaBetweenness"
                    name="Delta Betweenness"
                    domain={[-0.02, 0.25]}
                    tick={{ fill: '#A6AEB8', fontSize: 11, fontFamily: 'JetBrains Mono' }}
                    tickLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    axisLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    label={{ 
                      value: 'Delta Betweenness', 
                      position: 'bottom', 
                      offset: 25,
                      fill: '#A6AEB8',
                      fontSize: 12,
                      fontFamily: 'Inter'
                    }}
                  />
                  
                  <YAxis
                    type="number"
                    dataKey="deltaEigenvector"
                    name="Delta Eigenvector"
                    domain={[-0.01, 0.14]}
                    tick={{ fill: '#A6AEB8', fontSize: 11, fontFamily: 'JetBrains Mono' }}
                    tickLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    axisLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    label={{ 
                      value: 'Delta Eigenvector', 
                      angle: -90, 
                      position: 'insideLeft',
                      offset: 10,
                      fill: '#A6AEB8',
                      fontSize: 12,
                      fontFamily: 'Inter'
                    }}
                  />
                  
                  <ZAxis type="number" dataKey="ensembleScore" range={[30, 150]} />
                  
                  <Tooltip 
                    content={<CustomTooltip />}
                    cursor={{ stroke: 'rgba(0, 240, 255, 0.3)', strokeWidth: 1, strokeDasharray: '4 4' }}
                  />
                  
                  {/* Reference lines at zero */}
                  <ReferenceLine 
                    x={0} 
                    stroke="rgba(255, 255, 255, 0.15)" 
                    strokeDasharray="4 4"
                  />
                  <ReferenceLine 
                    y={0} 
                    stroke="rgba(255, 255, 255, 0.15)" 
                    strokeDasharray="4 4"
                  />
                  
                  {/* Other genes (small, gray) */}
                  <Scatter
                    name="Other"
                    data={otherGenes}
                    fill="rgba(166, 174, 184, 0.4)"
                    stroke="none"
                    onMouseEnter={(data) => setHoveredPoint(data)}
                    onMouseLeave={() => setHoveredPoint(null)}
                  />
                  
                  {/* Candidate genes (cyan) */}
                  <Scatter
                    name="Candidates"
                    data={candidateGenes}
                    fill="#00F0FF"
                    fillOpacity={0.7}
                    stroke="#00F0FF"
                    strokeWidth={1}
                    onMouseEnter={(data) => setHoveredPoint(data)}
                    onMouseLeave={() => setHoveredPoint(null)}
                  />
                  
                  {/* Top genes (gold, with glow) */}
                  <Scatter
                    name="Top"
                    data={topGenes}
                    fill="#FFD166"
                    fillOpacity={0.9}
                    stroke="#FFD166"
                    strokeWidth={2}
                    onMouseEnter={(data) => setHoveredPoint(data)}
                    onMouseLeave={() => setHoveredPoint(null)}
                  />
                </ScatterChart>
              </ResponsiveContainer>
            </div>

            {/* Quadrant Labels */}
            <div className="relative h-8 mt-2">
              <span className="absolute left-[25%] text-xs text-[#A6AEB8]/50">
                Decreased Influence
              </span>
              <span className="absolute right-[15%] text-xs text-[#00F0FF]/70">
                Increased Influence →
              </span>
            </div>
          </motion.div>
        </div>
      </div>
    </section>
  );
};

export default ScatterPlotSection;
