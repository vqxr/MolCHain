import React, { useState, useEffect } from 'react';
import Dashboard from './components/Dashboard';
import './App.css';

const API_URL = '/api';  // relative path; Vite will proxy to backend in dev

function App() {
  const [molecules, setMolecules] = useState([]);
  const [connected, setConnected] = useState(false);
  const [loading, setLoading] = useState(true);
  const [stats, setStats] = useState(null);
  const [aiAgents, setAiAgents] = useState([]);
  const [flarePrice, setFlarePrice] = useState(null);
  const [showDashboard, setShowDashboard] = useState(true);

  // ========== ALL FUNCTIONS DEFINED HERE ==========
  
  // Function 1: Test Flare Oracle
  const testFlare = () => {
    fetch(`${API_URL}/flare/aspirin_001`)
      .then(res => res.json())
      .then(data => {
        alert(`Flare Oracle Value: $${data.value_usd}\nConfidence: ${Math.round(data.confidence * 100)}%`);
      })
      .catch(err => {
        alert('Flare API error: ' + err.message);
      });
  };

  // Function 2: Run AI Analysis
  const runAIAnalysis = (moleculeId) => {
    const moleculeName = molecules.find(m => m.id === moleculeId)?.name || 'Unknown';
    
    fetch(`${API_URL}/analyze`, {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ molecule_id: moleculeId })
    })
    .then(res => res.json())
    .then(data => {
      if (data.success) {
        alert(`✅ AI Analysis Complete for ${moleculeName}!\nConsensus Score: ${data.consensus_score}%\nBlockchain TX: ${data.blockchain.tx_hash}`);
        // Refresh molecules
        fetch(`${API_URL}/molecules`)
          .then(res => res.json())
          .then(data => {
            if (data.success) setMolecules(data.molecules);
          });
      }
    })
    .catch(err => {
      alert('AI Analysis failed: ' + err.message);
    });
  };

  // Function 3: Simulate AI Analysis (general)
  const simulateAIAnalysis = () => {
    alert('🤖 AI Agents analyzing... This would call Python AI backend!');
  };

  // ========== USE EFFECT ==========
  useEffect(() => {
  const loadData = async () => {
    try {
      // First check if backend is alive
        const healthRes = await fetch(`${API_URL}/health`);
        const text = await healthRes.text();
        console.log('Raw response:', text);
        const healthData = JSON.parse(text);
        console.log('Parsed JSON:', healthData);
      
      if (healthData.status === 'healthy') {
        setConnected(true);
        
        // Load molecules (essential)
        try {
          const molRes = await fetch(`${API_URL}/molecules`);
          const molData = await molRes.json();
          if (molData.success) setMolecules(molData.molecules || []);
        } catch (err) {
          console.log('Could not load molecules:', err);
        }
        
        // Load stats (optional - if fails, continue)
        try {
          const statsRes = await fetch(`${API_URL}/stats`);
          const statsData = await statsRes.json();
          if (statsData.success) setStats(statsData.stats);
        } catch (err) {
          console.log('Could not load stats:', err);
        }
        
        // Load agents (optional)
        try {
          const agentsRes = await fetch(`${API_URL}/agents`);
          const agentsData = await agentsRes.json();
          if (agentsData.success) setAiAgents(agentsData.swarm?.agents || []);
        } catch (err) {
          console.log('Could not load agents:', err);
        }
        
        // Load flare price (optional)
        try {
          const priceRes = await fetch(`${API_URL}/flare/price`);
          const priceData = await priceRes.json();
          if (priceData.success) setFlarePrice(priceData);
        } catch (err) {
          console.log('Could not load flare price:', err);
          // Set a default price if API fails
          setFlarePrice({
            success: true,
            price_usd: 0.035,
            confidence: 0.95,
            oracle: "Flare FTSO"
          });
        }
      }
    } catch (err) {
      console.log('Cannot connect to backend:', err);
      setConnected(false);
    } finally {
      setLoading(false);
    }
  };

  loadData();
}, []);

  // If Dashboard is enabled, show it as the main view
  if (showDashboard) {
    return <Dashboard />;
  }

  return (
    <div style={{
      minHeight: '100vh',
      width: '100vw',  
      backgroundColor: '#0f172a',
      color: 'white',
      fontFamily: 'system-ui, -apple-system, sans-serif',
      padding: '20px',
      overflowX: 'hidden'  
    }}>
      {/* Header */}
      <div style={{
        maxWidth: '1200px',
        margin: '0 auto'
      }}>
        <div style={{
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          marginBottom: '40px'
        }}>
          <div>
            <h1 style={{
              fontSize: '36px',
              fontWeight: 'bold',
              marginBottom: '8px'
            }}>🧬 MolChain</h1>
            <p style={{
              color: '#94a3b8',
              fontSize: '18px'
            }}>AI Agent Swarms + Flare Blockchain</p>
          </div>
          
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: '12px'
          }}>
            <div style={{
              padding: '8px 16px',
              backgroundColor: connected ? '#10b98120' : '#f59e0b20',
              color: connected ? '#10b981' : '#f59e0b',
              borderRadius: '20px',
              fontSize: '14px',
              border: `1px solid ${connected ? '#10b98140' : '#f59e0b40'}`
            }}>
              {connected ? '✅ Backend Online' : '⚠️ Backend Offline'}
            </div>
            <button style={{
              padding: '10px 20px',
              backgroundColor: '#3b82f6',
              color: 'white',
              border: 'none',
              borderRadius: '8px',
              fontWeight: '500',
              cursor: 'pointer'
            }}>
              Connect Wallet
            </button>
          </div>
        </div>

        {/* Stats */}
        <div style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(4, 1fr)',
          gap: '16px',
          marginBottom: '40px'
        }}>
          {[
            { label: 'AI Agents', value: '5', emoji: '🤖' },
            { label: 'On-Chain', value: connected ? (stats?.verified_molecules || molecules.filter(m => m.on_blockchain).length) : '0', emoji: '🔗' },
            { label: 'Molecules', value: stats?.total_molecules || molecules.length, emoji: '🧪' },
            { label: 'Avg Score', value: stats?.avg_ai_score ? stats.avg_ai_score + '%' : (molecules.length > 0 ? Math.round(molecules.reduce((sum, m) => sum + (m.ai_score || 0), 0) / molecules.length) + '%' : '0%'), emoji: '📊' }
          ].map((stat, i) => (
            <div key={i} style={{
              backgroundColor: '#1e293b',
              padding: '20px',
              borderRadius: '12px',
              border: '1px solid #334155'
            }}>
              <div style={{
                fontSize: '24px',
                marginBottom: '8px'
              }}>{stat.emoji}</div>
              <div style={{
                fontSize: '28px',
                fontWeight: 'bold',
                marginBottom: '4px'
              }}>{stat.value}</div>
              <div style={{
                color: '#94a3b8',
                fontSize: '14px'
              }}>{stat.label}</div>
            </div>
          ))}
        </div>

        {/* Flare Price Card */}
        {flarePrice && (
          <div style={{
            backgroundColor: 'rgba(168, 85, 247, 0.1)',
            border: '1px solid rgba(168, 85, 247, 0.3)',
            borderRadius: '12px',
            padding: '20px',
            marginBottom: '24px',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between'
          }}>
            <div>
              <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
                <span style={{ fontSize: '24px' }}>🌐</span>
                <h3 style={{ fontSize: '18px', fontWeight: 'bold' }}>Flare FTSO Oracle</h3>
              </div>
              <p style={{ color: '#cbd5e1', fontSize: '14px' }}>Live price feed from blockchain</p>
            </div>
            <div style={{ textAlign: 'right' }}>
              <div style={{ fontSize: '32px', fontWeight: 'bold', color: '#a78bfa' }}>
                ${flarePrice.price_usd}
              </div>
              <div style={{ fontSize: '14px', color: '#cbd5e1' }}>
                Confidence: {Math.round(flarePrice.confidence * 100)}%
              </div>
            </div>
          </div>
        )}

        {/* AI Agents Card */}
        {aiAgents.length > 0 && (
          <div style={{
            backgroundColor: '#1e293b',
            border: '1px solid #334155',
            borderRadius: '12px',
            padding: '20px',
            marginBottom: '24px'
          }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px' }}>
              <h3 style={{ fontSize: '18px', fontWeight: 'bold', display: 'flex', alignItems: 'center', gap: '8px' }}>
                <span>🤖</span> AI Agent Swarm
              </h3>
              <div style={{
                padding: '4px 12px',
                backgroundColor: '#10b98120',
                color: '#10b981',
                borderRadius: '20px',
                fontSize: '12px',
                display: 'flex',
                alignItems: 'center',
                gap: '6px'
              }}>
                <div style={{ width: '6px', height: '6px', backgroundColor: '#10b981', borderRadius: '50%' }}></div>
                All Active
              </div>
            </div>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(5, 1fr)', gap: '12px' }}>
              {aiAgents.map(agent => (
                <div key={agent.name} style={{ textAlign: 'center' }}>
                  <div style={{ fontSize: '20px', marginBottom: '8px' }}>⚡</div>
                  <div style={{ fontSize: '12px', fontWeight: 'bold', whiteSpace: 'nowrap' }}>
                    {agent.name.split(' ')[0]}
                  </div>
                  <div style={{ fontSize: '11px', color: '#94a3b8' }}>
                    {Math.round(agent.accuracy * 100)}% acc
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Buttons */}
        <div style={{
          display: 'flex',
          gap: '12px',
          marginBottom: '40px'
        }}>
          <button 
            onClick={testFlare}
            style={{
              padding: '12px 24px',
              backgroundColor: '#8b5cf6',
              color: 'white',
              border: 'none',
              borderRadius: '8px',
              fontWeight: '500',
              cursor: 'pointer',
              fontSize: '16px'
            }}
          >
            ⚡ Test Flare Oracle
          </button>
          <button 
            onClick={simulateAIAnalysis}
            style={{
              padding: '12px 24px',
              backgroundColor: '#10b981',
              color: 'white',
              border: 'none',
              borderRadius: '8px',
              fontWeight: '500',
              cursor: 'pointer',
              fontSize: '16px'
            }}
          >
            🧠 Simulate AI Analysis
          </button>
        </div>

        {/* Molecules Section */}
        <div>
          <h2 style={{
            fontSize: '24px',
            fontWeight: 'bold',
            marginBottom: '24px'
          }}>Molecular Database</h2>
          
          {loading ? (
            <div style={{
              textAlign: 'center',
              padding: '60px 20px'
            }}>
              <div style={{
                width: '50px',
                height: '50px',
                border: '3px solid #3b82f6',
                borderTopColor: 'transparent',
                borderRadius: '50%',
                animation: 'spin 1s linear infinite',
                margin: '0 auto 20px'
              }}></div>
              <p>Connecting to Python backend...</p>
              <p style={{ color: '#64748b', marginTop: '8px' }}>
                Make sure backend is running on port 5000
              </p>
            </div>
          ) : molecules.length === 0 ? (
            <div style={{
              textAlign: 'center',
              padding: '60px 20px',
              backgroundColor: '#1e293b',
              borderRadius: '12px',
              border: '1px solid #334155'
            }}>
              <div style={{
                fontSize: '48px',
                marginBottom: '20px'
              }}>⚠️</div>
              <h3 style={{
                fontSize: '20px',
                marginBottom: '12px'
              }}>No Molecules Found</h3>
              <p style={{
                color: '#94a3b8',
                marginBottom: '24px',
                maxWidth: '500px',
                margin: '0 auto 24px'
              }}>
                {connected 
                  ? 'Backend is connected but returned no molecules. Check your Python API.'
                  : 'Cannot connect to Python backend. Make sure it\'s running.'}
              </p>
              <div style={{
                backgroundColor: '#0f172a',
                padding: '16px',
                borderRadius: '8px',
                fontFamily: 'monospace',
                fontSize: '14px',
                color: '#94a3b8',
                display: 'inline-block'
              }}>
                $ cd backend && python3 app.py
              </div>
            </div>
          ) : (
            <div style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fill, minmax(300px, 1fr))',
              gap: '20px'
            }}>
              {molecules.map(mol => (
                <div key={mol.id} style={{
                  backgroundColor: '#1e293b',
                  padding: '24px',
                  borderRadius: '12px',
                  border: '1px solid #334155'
                }}>
                  <div style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'flex-start',
                    marginBottom: '20px'
                  }}>
                    <div>
                      <h3 style={{
                        fontSize: '20px',
                        fontWeight: 'bold',
                        margin: '0 0 8px 0'
                      }}>{mol.name}</h3>
                      <p style={{
                        color: '#94a3b8',
                        fontFamily: 'monospace',
                        margin: '0'
                      }}>{mol.formula}</p>
                    </div>
                    <div style={{
                      padding: '6px 12px',
                      backgroundColor: mol.on_blockchain ? '#10b98120' : '#f59e0b20',
                      color: mol.on_blockchain ? '#10b981' : '#f59e0b',
                      borderRadius: '20px',
                      fontSize: '12px',
                      border: `1px solid ${mol.on_blockchain ? '#10b98140' : '#f59e0b40'}`
                    }}>
                      {mol.on_blockchain ? '🔗 On-Chain' : '⚡ Analyzing'}
                    </div>
                  </div>
                  
                  {/* Score bar */}
                  <div style={{ marginBottom: '20px' }}>
                    <div style={{
                      display: 'flex',
                      justifyContent: 'space-between',
                      marginBottom: '8px',
                      color: '#94a3b8'
                    }}>
                      <span>AI Consensus Score</span>
                      <span style={{ fontWeight: 'bold', color: 'white' }}>{mol.ai_score}%</span>
                    </div>
                    <div style={{
                      height: '8px',
                      backgroundColor: '#334155',
                      borderRadius: '4px',
                      overflow: 'hidden'
                    }}>
                      <div style={{
                        width: `${mol.ai_score}%`,
                        height: '100%',
                        background: 'linear-gradient(90deg, #3b82f6, #8b5cf6)',
                        borderRadius: '4px'
                      }}></div>
                    </div>
                  </div>
                  
                  {/* Buttons */}
                  <div style={{ marginTop: '16px' }}>
                    <button
                      onClick={() => runAIAnalysis(mol.id)}
                      style={{
                        width: '100%',
                        padding: '8px',
                        backgroundColor: 'rgba(59, 130, 246, 0.1)',
                        color: '#3b82f6',
                        border: '1px solid rgba(59, 130, 246, 0.3)',
                        borderRadius: '6px',
                        cursor: 'pointer',
                        fontSize: '12px',
                        fontWeight: '500',
                        marginBottom: '8px'
                      }}
                    >
                      🤖 Run AI Analysis
                    </button>
                    
                    <a
                      href={`/molecule/${mol.id}`}
                      style={{
                        display: 'block',
                        width: '100%',
                        padding: '8px',
                        backgroundColor: 'rgba(168, 85, 247, 0.1)',
                        color: '#a78bfa',
                        border: '1px solid rgba(168, 85, 247, 0.3)',
                        borderRadius: '6px',
                        textDecoration: 'none',
                        textAlign: 'center',
                        fontSize: '12px',
                        fontWeight: '500'
                      }}
                    >
                      🔍 View Details
                    </a>
                  </div>
                  
                  {/* Transaction hash */}
                  {mol.tx_hash && (
                    <div style={{
                      backgroundColor: '#0f172a',
                      padding: '12px',
                      borderRadius: '8px',
                      fontFamily: 'monospace',
                      fontSize: '12px',
                      color: '#94a3b8',
                      marginTop: '16px'
                    }}>
                      TX: {mol.tx_hash}
                    </div>
                  )}
                </div>
              ))}
            </div>
          )}
        </div>

        {/* Footer */}
        <div style={{
          marginTop: '60px',
          paddingTop: '40px',
          borderTop: '1px solid #334155',
          color: '#64748b',
          textAlign: 'center'
        }}>
          <p>🚀 Built for ETH Oxford 2026 Hackathon</p>
          <p style={{ fontSize: '14px', marginTop: '8px' }}>
            Python Flask Backend • React Frontend • Flare Blockchain • AI Agent Swarms
          </p>
        </div>
      </div>

      {/* Add spinner animation */}
      <style>{`
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
      `}</style>
    </div>
  );
}

export default App;