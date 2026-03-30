import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import NetworkGraph from '../components/NetworkGraph';
import { topRewiredGenes } from '../data/mockData';

const NetworkExplorer = () => {
  const [tissue, setTissue] = useState('normal');
  const [selectedGene, setSelectedGene] = useState(null);
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);

  const geneDetails = selectedGene 
    ? topRewiredGenes.find(g => g.gene === selectedGene)
    : null;

  const handleNodeClick = (geneId) => {
    setSelectedGene(geneId);
    setIsSidebarOpen(!!geneId);
  };

  const handleReset = () => {
    setSelectedGene(null);
    setIsSidebarOpen(false);
  };

  return (
    <section id="explorer" className="relative w-full min-h-screen bg-[#0A0A0A] py-12 px-4 sm:px-8">
      {/* Background gradient */}
      <div className="absolute inset-0 gradient-mesh pointer-events-none opacity-30" />

      <div className="relative z-10 max-w-[1800px] mx-auto">
        {/* Header */}
        <motion.div
          className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4 mb-8"
          initial={{ opacity: 0, y: -20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.5 }}
        >
          <div>
            <h2 className="font-display text-[clamp(24px,3vw,36px)] font-bold text-[#F2F4F8]">
              Network Explorer
            </h2>
            <p className="mt-1 text-sm text-[#A6AEB8]">
              Visualize gene interactions and identify rewiring patterns
            </p>
          </div>

          {/* Tissue Toggle */}
          <div className="toggle-pill">
            <div 
              className="toggle-slider"
              style={{
                left: tissue === 'normal' ? '4px' : 'calc(50% + 2px)',
                background: tissue === 'normal' 
                  ? 'linear-gradient(135deg, rgba(0, 240, 255, 0.3), rgba(0, 240, 255, 0.1))'
                  : 'linear-gradient(135deg, rgba(255, 45, 141, 0.3), rgba(255, 45, 141, 0.1))',
                border: `1px solid ${tissue === 'normal' ? 'rgba(0, 240, 255, 0.4)' : 'rgba(255, 45, 141, 0.4)'}`,
              }}
            />
            <button
              onClick={() => setTissue('normal')}
              className={`toggle-option py-2 ${tissue === 'normal' ? 'text-[#F2F4F8]' : 'text-[#A6AEB8]'}`}
            >
              <span className="flex items-center justify-center gap-2">
                <span 
                  className="w-2 h-2 rounded-full"
                  style={{ backgroundColor: '#00F0FF', boxShadow: '0 0 8px #00F0FF' }}
                />
                Normal
              </span>
            </button>
            <button
              onClick={() => setTissue('tumor')}
              className={`toggle-option py-2 ${tissue === 'tumor' ? 'text-[#F2F4F8]' : 'text-[#A6AEB8]'}`}
            >
              <span className="flex items-center justify-center gap-2">
                <span 
                  className="w-2 h-2 rounded-full"
                  style={{ backgroundColor: '#FF2D8D', boxShadow: '0 0 8px #FF2D8D' }}
                />
                Tumor
              </span>
            </button>
          </div>
        </motion.div>

        {/* Main Content Area */}
        <div className="flex flex-col lg:flex-row gap-6">
          {/* Network Graph */}
          <motion.div
            className="flex-1 glass-card network-container"
            style={{ height: '70vh', minHeight: '500px' }}
            initial={{ opacity: 0, x: -40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true }}
            transition={{ duration: 0.6, delay: 0.1 }}
          >
            <NetworkGraph
              tissue={tissue}
              selectedGene={selectedGene}
              onNodeClick={handleNodeClick}
              width="100%"
              height="100%"
            />
          </motion.div>

          {/* Sidebar */}
          <AnimatePresence mode="wait">
            {isSidebarOpen && geneDetails && (
              <motion.div
                className="w-full lg:w-80 glass-card p-6"
                initial={{ opacity: 0, x: 40 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, x: 40 }}
                transition={{ duration: 0.3, ease: [0.25, 0.46, 0.45, 0.94] }}
              >
                {/* Gene Header */}
                <div className="flex items-center justify-between mb-6">
                  <div>
                    <h3 className="font-display text-2xl font-bold text-[#FFD166]">
                      {geneDetails.gene}
                    </h3>
                    <span className="text-xs text-[#A6AEB8] uppercase tracking-wider">
                      Rank #{geneDetails.rank}
                    </span>
                  </div>
                  <button
                    onClick={handleReset}
                    className="w-8 h-8 rounded-lg border border-white/10 flex items-center justify-center
                               text-[#A6AEB8] hover:text-[#F2F4F8] hover:border-white/20 transition-colors"
                  >
                    <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                    </svg>
                  </button>
                </div>

                {/* Ensemble Score */}
                <div className="mb-6">
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-sm text-[#A6AEB8]">Ensemble Score</span>
                    <span className="font-mono-data text-lg font-bold text-[#00F0FF]">
                      {geneDetails.ensembleScore.toFixed(2)}
                    </span>
                  </div>
                  <div className="progress-bar">
                    <div 
                      className="progress-bar bg-gradient-to-r from-[#00F0FF] to-[#00F0FF]/50"
                      style={{ width: `${geneDetails.ensembleScore * 100}%` }}
                    />
                  </div>
                </div>

                {/* Centrality Scores */}
                <div className="space-y-4 mb-6">
                  <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider">
                    Centrality Scores
                  </h4>
                  
                  {/* Betweenness */}
                  <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs text-[#A6AEB8]">Betweenness (Delta)</span>
                      <span 
                        className="font-mono-data text-sm font-bold"
                        style={{ color: geneDetails.deltaBetweenness > 0 ? '#00F0FF' : '#FF2D8D' }}
                      >
                        {geneDetails.deltaBetweenness > 0 ? '+' : ''}
                        {geneDetails.deltaBetweenness.toFixed(2)}
                      </span>
                    </div>
                    <div className="flex items-center gap-2 text-xs text-[#A6AEB8]">
                      <span>N: {geneDetails.normalBetweenness.toFixed(2)}</span>
                      <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 8l4 4m0 0l-4 4m4-4H3" />
                      </svg>
                      <span>T: {geneDetails.tumorBetweenness.toFixed(2)}</span>
                    </div>
                  </div>

                  {/* Eigenvector */}
                  <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs text-[#A6AEB8]">Eigenvector (Delta)</span>
                      <span 
                        className="font-mono-data text-sm font-bold"
                        style={{ color: geneDetails.deltaEigenvector > 0 ? '#00F0FF' : '#FF2D8D' }}
                      >
                        {geneDetails.deltaEigenvector > 0 ? '+' : ''}
                        {geneDetails.deltaEigenvector.toFixed(2)}
                      </span>
                    </div>
                    <div className="flex items-center gap-2 text-xs text-[#A6AEB8]">
                      <span>N: {geneDetails.normalEigenvector.toFixed(2)}</span>
                      <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 8l4 4m0 0l-4 4m4-4H3" />
                      </svg>
                      <span>T: {geneDetails.tumorEigenvector.toFixed(2)}</span>
                    </div>
                  </div>
                </div>

                {/* Description */}
                <div className="mb-6">
                  <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-2">
                    Description
                  </h4>
                  <p className="text-sm text-[#A6AEB8] leading-relaxed">
                    {geneDetails.description}
                  </p>
                </div>

                {/* Actions */}
                <div className="flex flex-col gap-2">
                  <button className="w-full py-3 px-4 rounded-xl bg-[#00F0FF]/10 border border-[#00F0FF]/30
                                     text-[#00F0FF] text-sm font-medium hover:bg-[#00F0FF]/20 transition-colors
                                     flex items-center justify-center gap-2">
                    <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 20l4-16m2 16l4-16M6 9h14M4 15h14" />
                    </svg>
                    Highlight Neighbors
                  </button>
                  <button className="w-full py-3 px-4 rounded-xl bg-white/5 border border-white/10
                                     text-[#F2F4F8] text-sm font-medium hover:bg-white/10 transition-colors
                                     flex items-center justify-center gap-2">
                    <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
                    </svg>
                    Export Subnetwork
                  </button>
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Placeholder when no gene selected */}
          {!isSidebarOpen && (
            <motion.div
              className="hidden lg:flex w-80 glass-card p-6 items-center justify-center"
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              transition={{ duration: 0.3 }}
            >
              <div className="text-center">
                <div className="w-16 h-16 mx-auto mb-4 rounded-full bg-white/[0.03] border border-white/[0.08]
                                flex items-center justify-center">
                  <svg className="w-8 h-8 text-[#A6AEB8]/50" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
                          d="M13 10V3L4 14h7v7l9-11h-7z" />
                  </svg>
                </div>
                <p className="text-sm text-[#A6AEB8]">
                  Click on a node in the network to view gene details
                </p>
              </div>
            </motion.div>
          )}
        </div>
      </div>
    </section>
  );
};

export default NetworkExplorer;
