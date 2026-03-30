import { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import NetworkGraph from '../components/NetworkGraph';

const HubRewiring = () => {
  const [geneInput, setGeneInput] = useState('');
  const [selectedGene, setSelectedGene] = useState(null);
  const [tissue, setTissue] = useState('normal');
  const [networkData, setNetworkData] = useState(null);
  const [initialNetwork, setInitialNetwork] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  // Load initial skeleton network on mount
  useEffect(() => {
    fetchInitialNetwork();
  }, []);

  const fetchInitialNetwork = async () => {
    try {
      const response = await fetch('http://localhost:8000/api/network/initial');
      if (!response.ok) throw new Error('Failed to fetch initial network');
      const data = await response.json();
      setInitialNetwork(data);
      setNetworkData(data);
    } catch (err) {
      console.error('Error loading initial network:', err);
      setError(err.message);
    }
  };

  const handleAnalyze = async () => {
    if (!geneInput.trim()) return;

    setLoading(true);
    setError(null);

    try {
      const response = await fetch(`http://localhost:8000/api/rewiring?gene=${geneInput.trim()}`);
      if (!response.ok) {
        throw new Error('Failed to fetch rewiring data');
      }
      const data = await response.json();
      setNetworkData(data);
      setSelectedGene(geneInput.trim());
    } catch (err) {
      setError(err.message);
      setNetworkData(initialNetwork); // Fall back to initial network
    } finally {
      setLoading(false);
    }
  };

  const handleNodeClick = (geneId) => {
    if (geneId) {
      setGeneInput(geneId);
      setSelectedGene(geneId);
      // Auto-fetch rewiring data for clicked gene
      fetchRewiringForGene(geneId);
    }
  };

  const fetchRewiringForGene = async (gene) => {
    setLoading(true);
    setError(null);

    try {
      const response = await fetch(`http://localhost:8000/api/rewiring?gene=${gene}`);
      if (!response.ok) throw new Error('Failed to fetch rewiring data');
      const data = await response.json();
      setNetworkData(data);
    } catch (err) {
      setError(err.message);
      setNetworkData(initialNetwork);
    } finally {
      setLoading(false);
    }
  };

  const calculateStats = () => {
    if (!networkData?.links) return { total: 0, gain: 0, loss: 0, maintained: 0 };

    const stats = networkData.links.reduce((acc, link) => {
      acc.total++;
      if (link.status === 'gained') acc.gain++;
      else if (link.status === 'lost') acc.loss++;
      else if (link.status === 'maintained') acc.maintained++;
      return acc;
    }, { total: 0, gain: 0, loss: 0, maintained: 0 });

    return stats;
  };

  const stats = calculateStats();
  const mode = selectedGene && networkData?.nodes?.some(n => n.type === 'hub') ? 'ego' : 'full';

  return (
    <div className="h-full flex flex-col p-6 gap-6">
      {/* Controls Panel */}
      <div className="flex-none">
        <div className="glass-card p-6">
          <div className="flex items-end gap-4">
            <div className="flex-1">
              <label className="block text-sm text-[#A6AEB8] mb-2">Gene Symbol</label>
              <input
                type="text"
                value={geneInput}
                onChange={(e) => setGeneInput(e.target.value.toUpperCase())}
                onKeyPress={(e) => e.key === 'Enter' && handleAnalyze()}
                placeholder="e.g., PARK2, TP53, EGFR"
                className="w-full px-4 py-3 rounded-xl bg-white/[0.03] border border-white/[0.08]
                           text-[#F2F4F8] placeholder-[#A6AEB8]/50
                           focus:outline-none focus:border-[#00F0FF]/50 focus:ring-1 focus:ring-[#00F0FF]/30
                           transition-all duration-200"
              />
            </div>

            <button
              onClick={handleAnalyze}
              disabled={loading || !geneInput.trim()}
              className="px-6 py-3 rounded-xl bg-[#00F0FF]/10 border border-[#00F0FF]/40
                         text-[#00F0FF] font-medium hover:bg-[#00F0FF]/20 hover:border-[#00F0FF]/60
                         disabled:opacity-50 disabled:cursor-not-allowed
                         transition-all duration-200 flex items-center gap-2"
            >
              {loading ? (
                <>
                  <svg className="w-5 h-5 animate-spin" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
                  </svg>
                  Analyzing...
                </>
              ) : (
                <>
                  <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                  </svg>
                  Analyze Hub
                </>
              )}
            </button>

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
          </div>

          {error && (
            <div className="mt-4 p-3 rounded-lg bg-[#FF2D8D]/10 border border-[#FF2D8D]/30 text-[#FF2D8D] text-sm">
              {error}
            </div>
          )}
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 flex gap-6 min-h-0">
        {/* Network Graph */}
        <div className="flex-1 glass-card network-container">
          <NetworkGraph
            tissue={tissue}
            selectedGene={selectedGene}
            onNodeClick={handleNodeClick}
            width="100%"
            height="100%"
            networkData={networkData}
            mode={mode}
          />
        </div>

        {/* Stats Panel */}
        {selectedGene && mode === 'ego' && (
          <motion.div
            className="w-80 glass-card p-6 flex flex-col gap-4"
            initial={{ opacity: 0, x: 40 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.3 }}
          >
            <div>
              <h3 className="font-display text-2xl font-bold text-[#FFD166] mb-1">
                {selectedGene}
              </h3>
              <p className="text-xs text-[#A6AEB8] uppercase tracking-wider">
                Hub Rewiring Analysis
              </p>
            </div>

            <div className="space-y-3">
              <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                <div className="text-xs text-[#A6AEB8] uppercase tracking-wider mb-1">
                  Total Connections
                </div>
                <div className="font-mono-data text-3xl font-bold text-[#00F0FF]">
                  {stats.total}
                </div>
              </div>

              <div className="p-4 rounded-xl bg-white/[0.03] border border-white/[0.06]">
                <div className="text-xs text-[#A6AEB8] uppercase tracking-wider mb-1">
                  Maintained Edges
                </div>
                <div className="font-mono-data text-2xl font-bold text-[#9B59D6]">
                  {stats.maintained}
                </div>
                <div className="text-xs text-[#A6AEB8] mt-1">stable connections</div>
              </div>
            </div>

            <div className="pt-4 border-t border-white/10">
              <h4 className="text-sm font-medium text-[#F2F4F8] uppercase tracking-wider mb-2">
                Neighbors
              </h4>
              <div className="space-y-2 max-h-40 overflow-y-auto">
                {networkData?.nodes
                  ?.filter(n => n.type !== 'hub')
                  ?.slice(0, 8)
                  .map((node) => (
                    <div
                      key={node.id}
                      className="flex items-center justify-between p-2 rounded-lg bg-white/[0.02] hover:bg-white/[0.05] cursor-pointer transition-colors"
                      onClick={() => {
                        setGeneInput(node.id);
                        fetchRewiringForGene(node.id);
                      }}
                    >
                      <div className="flex items-center gap-2">
                        <div 
                          className="w-2 h-2 rounded-full"
                          style={{
                            backgroundColor: node.status === 'lost' ? '#00F0FF' : 
                                           node.status === 'gained' ? '#FF2D8D' : '#9B59D6'
                          }}
                        />
                        <span className="text-sm text-[#F2F4F8] font-mono-data">{node.id}</span>
                      </div>
                      <svg className="w-4 h-4 text-[#A6AEB8]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                      </svg>
                    </div>
                  ))}
              </div>
            </div>
          </motion.div>
        )}
      </div>
    </div>
  );
};

export default HubRewiring;