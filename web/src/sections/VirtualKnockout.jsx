import { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Cell,
  ZAxis
} from 'recharts';

const VirtualKnockout = () => {
  const [knockoutData, setKnockoutData] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    fetchKnockoutData();
  }, []);

  const fetchKnockoutData = async () => {
    setLoading(true);
    setError(null);

    try {
      const response = await fetch('http://localhost:8000/api/knockout');
      if (!response.ok) {
        throw new Error('Failed to fetch knockout data');
      }
      const data = await response.json();
      setKnockoutData(data);
    } catch (err) {
      setError(err.message);
      setKnockoutData([]);
    } finally {
      setLoading(false);
    }
  };

  const getBarColor = (vitality) => {
    if (vitality > 20) return '#FF2D8D';
    if (vitality > 10) return '#FFD166';
    return '#00F0FF';
  };

  const CustomTooltip = ({ active, payload }) => {
    if (active && payload && payload.length) {
      const data = payload[0].payload;
      return (
        <div className="glass-card-sm p-3 border border-[#00F0FF]/30">
          <p className="font-mono-data font-bold text-[#F2F4F8] mb-1 text-sm">{data.gene}</p>
          <div className="space-y-0.5 text-xs">
            <p className="text-[#A6AEB8]">
              Vitality: <span className="text-[#00F0FF] font-bold">{data.vitality_score.toFixed(2)}</span>
            </p>
            <p className="text-[#A6AEB8]">
              ASPL: <span className="text-[#FFD166]">{data.aspl_knockout.toFixed(1)}</span>
            </p>
          </div>
        </div>
      );
    }
    return null;
  };

  // Limit data for visualization - show top 50
  const displayData = knockoutData.slice(0, 50);
  
  // Prepare scatter plot data - no scaling needed as ASPL is already in 1700-1800 range
  const scatterData = displayData.map((item, idx) => ({
    ...item,
    rank: idx + 1,
    size: Math.max(30, 100 - idx * 1.5) // Size decreases with rank
  }));

  return (
    <div className="h-full flex flex-col p-6 gap-6 overflow-hidden">
      {/* Header */}
      <div className="flex-none">
        <div className="glass-card p-6">
          <div className="flex items-center justify-between">
            <div>
              <h2 className="font-display text-2xl font-bold text-[#F2F4F8] mb-1">
                Virtual Knockout Analysis
              </h2>
              <p className="text-sm text-[#A6AEB8]">
                Hub genes critical for network integrity. Higher vitality = more essential.
              </p>
            </div>

            <button
              onClick={fetchKnockoutData}
              disabled={loading}
              className="px-4 py-2 rounded-xl bg-white/5 border border-white/10
                         text-[#A6AEB8] hover:text-[#00F0FF] hover:border-[#00F0FF]/30
                         disabled:opacity-50 disabled:cursor-not-allowed
                         transition-all duration-200 flex items-center gap-2"
            >
              <svg className={`w-4 h-4 ${loading ? 'animate-spin' : ''}`} fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
              </svg>
              Refresh
            </button>
          </div>

          {error && (
            <div className="mt-4 p-3 rounded-lg bg-[#FF2D8D]/10 border border-[#FF2D8D]/30 text-[#FF2D8D] text-sm">
              {error}
            </div>
          )}
        </div>
      </div>

      {/* Main Content */}
      {loading ? (
        <div className="flex-1 flex items-center justify-center">
          <div className="text-center">
            <svg className="w-12 h-12 mx-auto mb-4 animate-spin text-[#00F0FF]" fill="none" viewBox="0 0 24 24">
              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
            </svg>
            <p className="text-sm text-[#A6AEB8]">Loading knockout analysis...</p>
          </div>
        </div>
      ) : knockoutData.length > 0 ? (
        <div className="flex-1 grid grid-cols-1 lg:grid-cols-2 gap-6 min-h-0 overflow-hidden">
          {/* Scatter Plot */}
          <motion.div
            className="glass-card p-6 flex flex-col overflow-hidden"
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.4 }}
          >
            <h3 className="font-display text-lg font-bold text-[#F2F4F8] mb-4 flex-none">
              Vitality vs Network Impact
            </h3>
            <div className="flex-1 min-h-0">
              <ResponsiveContainer width="100%" height="100%">
                <ScatterChart margin={{ top: 10, right: 20, bottom: 40, left: 50 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255, 255, 255, 0.06)" />
                  <XAxis
                    dataKey="aspl_knockout"
                    type="number"
                    name="ASPL After Knockout"
                    domain={['auto', 'auto']}
                    tick={{ fill: '#A6AEB8', fontSize: 10, fontFamily: 'JetBrains Mono' }}
                    tickLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    axisLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    label={{ 
                      value: 'ASPL After Knockout', 
                      position: 'bottom', 
                      offset: 20,
                      fill: '#A6AEB8',
                      fontSize: 11
                    }}
                  />
                  <YAxis
                    dataKey="vitality_score"
                    type="number"
                    name="Vitality Score"
                    domain={[0, 'auto']}
                    tick={{ fill: '#A6AEB8', fontSize: 10, fontFamily: 'JetBrains Mono' }}
                    tickLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    axisLine={{ stroke: 'rgba(255, 255, 255, 0.1)' }}
                    label={{ 
                      value: 'Vitality Score', 
                      angle: -90, 
                      position: 'insideLeft',
                      fill: '#A6AEB8',
                      fontSize: 11
                    }}
                  />
                  <ZAxis dataKey="size" range={[30, 150]} />
                  <Tooltip content={<CustomTooltip />} cursor={{ strokeDasharray: '3 3' }} />
                  <Scatter data={scatterData} shape="circle">
                    {scatterData.map((entry, index) => (
                      <Cell 
                        key={`cell-${index}`} 
                        fill={getBarColor(entry.vitality_score)}
                        fillOpacity={0.7}
                      />
                    ))}
                  </Scatter>
                </ScatterChart>
              </ResponsiveContainer>
            </div>
            <div className="flex-none mt-3 flex items-center justify-center gap-4 text-xs">
              <div className="flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-[#FF2D8D]" />
                <span className="text-[#A6AEB8]">Critical (&gt;20)</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-[#FFD166]" />
                <span className="text-[#A6AEB8]">High (10-20)</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-[#00F0FF]" />
                <span className="text-[#A6AEB8]">Moderate (&lt;10)</span>
              </div>
            </div>
          </motion.div>

          {/* Data Table with Scrolling */}
          <motion.div
            className="glass-card flex flex-col overflow-hidden"
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.4, delay: 0.1 }}
          >
            <div className="flex-none p-6 border-b border-white/10">
              <div className="flex items-center justify-between">
                <h3 className="font-display text-lg font-bold text-[#F2F4F8]">
                  Top Critical Genes
                </h3>
                <span className="text-xs text-[#A6AEB8]">
                  Showing {displayData.length} of {knockoutData.length}
                </span>
              </div>
            </div>
            <div className="flex-1 overflow-y-auto">
              <table className="w-full">
                <thead className="sticky top-0 bg-[#0A0A0A] z-10 border-b border-white/10">
                  <tr>
                    <th className="px-4 py-3 text-left text-xs font-medium text-[#A6AEB8] uppercase tracking-wider">
                      Rank
                    </th>
                    <th className="px-4 py-3 text-left text-xs font-medium text-[#A6AEB8] uppercase tracking-wider">
                      Gene
                    </th>
                    <th className="px-4 py-3 text-right text-xs font-medium text-[#A6AEB8] uppercase tracking-wider">
                      Vitality
                    </th>
                    <th className="px-4 py-3 text-right text-xs font-medium text-[#A6AEB8] uppercase tracking-wider">
                      ASPL Δ
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-white/5">
                  {displayData.map((gene, index) => (
                    <motion.tr
                      key={gene.gene}
                      className="hover:bg-white/[0.03] transition-colors"
                      initial={{ opacity: 0, y: 10 }}
                      animate={{ opacity: 1, y: 0 }}
                      transition={{ duration: 0.2, delay: index * 0.01 }}
                    >
                      <td className="px-4 py-3 whitespace-nowrap">
                        <span
                          className={`inline-flex items-center justify-center w-7 h-7 rounded-full text-xs font-bold ${
                            index < 3
                              ? 'bg-[#FF2D8D]/20 text-[#FF2D8D] border border-[#FF2D8D]/40'
                              : index < 10
                              ? 'bg-[#FFD166]/20 text-[#FFD166] border border-[#FFD166]/40'
                              : 'bg-white/5 text-[#A6AEB8]'
                          }`}
                        >
                          {index + 1}
                        </span>
                      </td>
                      <td className="px-4 py-3 whitespace-nowrap">
                        <span className="font-mono-data font-bold text-[#F2F4F8] text-sm">
                          {gene.gene}
                        </span>
                      </td>
                      <td className="px-4 py-3 whitespace-nowrap text-right">
                        <span
                          className="font-mono-data text-sm font-bold"
                          style={{ color: getBarColor(gene.vitality_score) }}
                        >
                          {gene.vitality_score.toFixed(2)}
                        </span>
                      </td>
                      <td className="px-4 py-3 whitespace-nowrap text-right">
                        <span className="font-mono-data text-sm text-[#A6AEB8]">
                          {gene.aspl_change.toFixed(2)}
                        </span>
                      </td>
                    </motion.tr>
                  ))}
                </tbody>
              </table>
            </div>
          </motion.div>
        </div>
      ) : (
        <div className="flex-1 flex items-center justify-center">
          <div className="text-center">
            <div className="w-20 h-20 mx-auto mb-6 rounded-full bg-white/[0.03] border border-white/[0.08] flex items-center justify-center">
              <svg className="w-10 h-10 text-[#A6AEB8]/50" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                      d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
              </svg>
            </div>
            <h3 className="text-lg font-medium text-[#F2F4F8] mb-2">No knockout data available</h3>
            <p className="text-sm text-[#A6AEB8]">Click refresh to load the analysis</p>
          </div>
        </div>
      )}
    </div>
  );
};

export default VirtualKnockout;