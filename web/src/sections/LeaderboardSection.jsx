import { useState } from 'react';
import { motion } from 'framer-motion';
import { topRewiredGenes } from '../data/mockData';

const ProgressBar = ({ value, color = '#00F0FF' }) => (
  <div className="progress-bar w-24">
    <motion.div
      className="progress-fill"
      style={{ backgroundColor: color }}
      initial={{ width: 0 }}
      whileInView={{ width: `${value * 100}%` }}
      viewport={{ once: true }}
      transition={{ duration: 0.8, delay: 0.3, ease: [0.25, 0.46, 0.45, 0.94] }}
    />
  </div>
);

const LeaderboardSection = ({ onGeneHighlight }) => {
  const [hoveredGene, setHoveredGene] = useState(null);

  const getRankColor = (rank) => {
    if (rank === 1) return '#FFD166';
    if (rank <= 3) return '#00F0FF';
    return '#A6AEB8';
  };

  const getRankIcon = (rank) => {
    if (rank === 1) return '🥇';
    if (rank === 2) return '🥈';
    if (rank === 3) return '🥉';
    return rank;
  };

  return (
    <section id="leaderboard" className="relative w-full min-h-screen bg-[#0A0A0A] py-20 px-4 sm:px-8 lg:px-16">
      {/* Background gradient */}
      <div className="absolute inset-0 gradient-mesh pointer-events-none opacity-30" />

      <div className="relative z-10 max-w-7xl mx-auto">
        {/* Header */}
        <motion.div
          className="mb-12"
          initial={{ opacity: 0, y: 30 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true, margin: "-100px" }}
          transition={{ duration: 0.6 }}
        >
          <h2 className="font-display text-[clamp(28px,3.5vw,44px)] font-bold text-[#F2F4F8]">
            Top rewiring candidates
          </h2>
          <p className="mt-3 text-[#A6AEB8] max-w-xl">
            Ranked by ensemble score. Hover a row to highlight the gene in the network, 
            or click to view detailed analysis.
          </p>
        </motion.div>

        {/* Leaderboard Table */}
        <motion.div
          className="glass-card overflow-hidden"
          initial={{ opacity: 0, y: 40 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true, margin: "-50px" }}
          transition={{ duration: 0.7, delay: 0.1 }}
        >
          <div className="overflow-x-auto">
            <table className="data-table w-full">
              <thead>
                <tr className="bg-white/[0.03]">
                  <th className="w-16 text-center">Rank</th>
                  <th>Gene</th>
                  <th>Ensemble Score</th>
                  <th>Delta Betweenness</th>
                  <th>Delta Eigenvector</th>
                  <th className="text-right">Action</th>
                </tr>
              </thead>
              <tbody>
                {topRewiredGenes.map((gene, index) => (
                  <motion.tr
                    key={gene.gene}
                    className="group cursor-pointer transition-colors duration-200"
                    style={{
                      backgroundColor: hoveredGene === gene.gene 
                        ? 'rgba(0, 240, 255, 0.06)' 
                        : 'transparent'
                    }}
                    initial={{ opacity: 0, x: -20 }}
                    whileInView={{ opacity: 1, x: 0 }}
                    viewport={{ once: true }}
                    transition={{ 
                      duration: 0.4, 
                      delay: 0.1 + index * 0.04,
                      ease: [0.25, 0.46, 0.45, 0.94]
                    }}
                    onMouseEnter={() => {
                      setHoveredGene(gene.gene);
                      onGeneHighlight?.(gene.gene);
                    }}
                    onMouseLeave={() => {
                      setHoveredGene(null);
                      onGeneHighlight?.(null);
                    }}
                  >
                    {/* Rank */}
                    <td className="text-center">
                      <span 
                        className="inline-flex items-center justify-center w-8 h-8 rounded-full text-sm font-bold"
                        style={{
                          backgroundColor: `${getRankColor(gene.rank)}20`,
                          color: getRankColor(gene.rank),
                          border: `1px solid ${getRankColor(gene.rank)}40`
                        }}
                      >
                        {getRankIcon(gene.rank)}
                      </span>
                    </td>

                    {/* Gene Name */}
                    <td>
                      <div className="flex items-center gap-3">
                        <span className="font-mono-data text-base font-bold text-[#F2F4F8]">
                          {gene.gene}
                        </span>
                        {gene.rank <= 3 && (
                          <span 
                            className="px-2 py-0.5 rounded text-[10px] font-medium uppercase tracking-wider"
                            style={{
                              backgroundColor: gene.rank === 1 
                                ? 'rgba(255, 209, 102, 0.15)' 
                                : 'rgba(0, 240, 255, 0.15)',
                              color: gene.rank === 1 ? '#FFD166' : '#00F0FF'
                            }}
                          >
                            {gene.rank === 1 ? 'Top Candidate' : 'High Priority'}
                          </span>
                        )}
                      </div>
                    </td>

                    {/* Ensemble Score */}
                    <td>
                      <div className="flex items-center gap-4">
                        <span className="font-mono-data text-sm font-bold text-[#00F0FF] w-12">
                          {gene.ensembleScore.toFixed(2)}
                        </span>
                        <ProgressBar 
                          value={gene.ensembleScore} 
                          color={gene.rank === 1 ? '#FFD166' : '#00F0FF'}
                        />
                      </div>
                    </td>

                    {/* Delta Betweenness */}
                    <td>
                      <span 
                        className="font-mono-data text-sm"
                        style={{ 
                          color: gene.deltaBetweenness > 0 ? '#00F0FF' : '#FF2D8D' 
                        }}
                      >
                        {gene.deltaBetweenness > 0 ? '+' : ''}
                        {gene.deltaBetweenness.toFixed(2)}
                      </span>
                    </td>

                    {/* Delta Eigenvector */}
                    <td>
                      <span 
                        className="font-mono-data text-sm"
                        style={{ 
                          color: gene.deltaEigenvector > 0 ? '#00F0FF' : '#FF2D8D' 
                        }}
                      >
                        {gene.deltaEigenvector > 0 ? '+' : ''}
                        {gene.deltaEigenvector.toFixed(2)}
                      </span>
                    </td>

                    {/* Action */}
                    <td className="text-right">
                      <button 
                        className="inline-flex items-center gap-2 px-4 py-2 rounded-lg text-sm
                                   bg-white/[0.03] border border-white/[0.08] text-[#A6AEB8]
                                   hover:text-[#00F0FF] hover:border-[#00F0FF]/30 hover:bg-[#00F0FF]/5
                                   transition-all duration-200 opacity-0 group-hover:opacity-100"
                        onClick={() => {
                          // Scroll to network explorer and select gene
                          document.getElementById('explorer')?.scrollIntoView({ behavior: 'smooth' });
                          onGeneHighlight?.(gene.gene);
                        }}
                      >
                        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} 
                                d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} 
                                d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
                        </svg>
                        Locate
                      </button>
                    </td>
                  </motion.tr>
                ))}
              </tbody>
            </table>
          </div>
        </motion.div>

        {/* Summary Cards */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mt-8">
          {[
            { 
              label: 'Analysis Method', 
              value: 'Ensemble Scoring',
              desc: 'Combines multiple centrality metrics'
            },
            { 
              label: 'Outlier Detection', 
              value: 'IQR-based',
              desc: 'Robust statistical filtering'
            },
            { 
              label: 'Export Format', 
              value: 'TSV/CSV/JSON',
              desc: 'Compatible with downstream tools'
            }
          ].map((item, i) => (
            <motion.div
              key={item.label}
              className="glass-card-sm p-5"
              initial={{ opacity: 0, y: 20 }}
              whileInView={{ opacity: 1, y: 0 }}
              viewport={{ once: true }}
              transition={{ duration: 0.5, delay: 0.5 + i * 0.1 }}
            >
              <div className="text-xs text-[#A6AEB8] uppercase tracking-wider mb-1">
                {item.label}
              </div>
              <div className="font-mono-data text-lg font-bold text-[#F2F4F8] mb-1">
                {item.value}
              </div>
              <div className="text-xs text-[#A6AEB8]/70">
                {item.desc}
              </div>
            </motion.div>
          ))}
        </div>
      </div>
    </section>
  );
};

export default LeaderboardSection;
