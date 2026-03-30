import MainApp from './sections/MainApp';

function App() {
  return (
    <div className="relative h-screen bg-[#0A0A0A] grain-overlay overflow-hidden">
      {/* Global Background Elements */}
      <div className="fixed inset-0 pointer-events-none z-0">
        <div
          className="absolute inset-0 opacity-30"
          style={{
            background: `
              radial-gradient(ellipse at 20% 20%, rgba(0, 240, 255, 0.05) 0%, transparent 50%),
              radial-gradient(ellipse at 80% 80%, rgba(255, 45, 141, 0.03) 0%, transparent 50%)
            `
          }}
        />
      </div>

      {/* Main App */}
      <MainApp />
    </div>
  );
}

export default App;
