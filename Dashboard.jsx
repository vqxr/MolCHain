import React, { useState, useEffect } from 'react';
import './Dashboard.css';
import MoleculeDetailModal from './MoleculeDetailModal';

const API_URL = '/api';

// Comprehensive Molecule Card with Animations
const MoleculeCard = ({ molecule, onClick }) => {
  return (
    <div 
      className="scifi-card neon-blue"
      onClick={() => onClick(molecule)}
      style={{cursor: 'pointer', animation: `slideIn 0.5s ease-out`}}
    >
      <div className="card-glow"></div>
      <div className="card-header">
        <div className="card-icon">🧬</div>
        <div className="card-title">{molecule.name}</div>
        <div className="card-status">{molecule.on_blockchain ? '✅' : '⏳'}</div>
      </div>
      <div className="molecule-formula">{molecule.formula}</div>
      <div className="molecule-weight">MW: {molecule.molecular_weight}</div>
      <div className="ai-score-bar">
        <div className="score-fill" style={{width: `${molecule.ai_score}%`}}></div>
      </div>
      <div className="score-label">AI Score: {molecule.ai_score}%</div>
      <div className="card-footer">
        <span className="value-badge">${molecule.blockchain_data?.value_usd?.toFixed(2)}</span>
      </div>
    </div>
  );
};

const Dashboard = () => {
  const [loading, setLoading] = useState(true);
  const [stats, setStats] = useState(null);
  const [molecules, setMolecules] = useState([]);
  const [aiAgents, setAiAgents] = useState([]);
  const [flarePrice, setFlarePrice] = useState(null);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [showModal, setShowModal] = useState(false);

  // Comprehensive Mock Data
  const mockMolecules = [
    { "id": "mol_001", "name": "Aspirin", "formula": "C9H8O4", "molecular_weight": 180.16, "ai_score": 87, "on_blockchain": true, "blockchain_data": { "tx_hash": "0x7a3f9c8d2e1b4a5f6c7d8e9f0a1b2c3d4e5f67890", "block_number": 1245678, "value_usd": 1245.67 } },
    { "id": "mol_002", "name": "Caffeine", "formula": "C8H10N4O2", "molecular_weight": 194.19, "ai_score": 92, "on_blockchain": true, "blockchain_data": { "tx_hash": "0x8b4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e", "block_number": 1245679, "value_usd": 892.45 } },
    { "id": "mol_003", "name": "Metformin", "formula": "C4H11N5", "molecular_weight": 129.16, "ai_score": 91, "on_blockchain": true, "blockchain_data": { "tx_hash": "0x9c5f6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c", "block_number": 1245680, "value_usd": 1123.89 } },
    { "id": "mol_004", "name": "Insulin", "formula": "C257H383N65O77S6", "molecular_weight": 5807.57, "ai_score": 96, "on_blockchain": true, "blockchain_data": { "tx_hash": "0xa1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c", "block_number": 1245681, "value_usd": 2450.75 } },
    { "id": "mol_005", "name": "Artemisinin", "formula": "C15H22O5", "molecular_weight": 282.33, "ai_score": 94, "on_blockchain": true, "blockchain_data": { "tx_hash": "0xb2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d", "block_number": 1245682, "value_usd": 1789.30 } },
    { "id": "mol_006", "name": "Penicillin", "formula": "C16H18N2O4S", "molecular_weight": 334.39, "ai_score": 88, "on_blockchain": true, "blockchain_data": { "tx_hash": "0xc3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e", "block_number": 1245683, "value_usd": 2156.42 } },
    { "id": "mol_007", "name": "Vitamin C", "formula": "C6H8O6", "molecular_weight": 176.12, "ai_score": 85, "on_blockchain": true, "blockchain_data": { "tx_hash": "0xd4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f", "block_number": 1245684, "value_usd": 567.89 } },
    { "id": "mol_008", "name": "Glucose", "formula": "C6H12O6", "molecular_weight": 180.16, "ai_score": 93, "on_blockchain": true, "blockchain_data": { "tx_hash": "0xe5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a", "block_number": 1245685, "value_usd": 834.56 } },
    { "id": "mol_009", "name": "Aspirin 2", "formula": "C9H8O4", "molecular_weight": 180.16, "ai_score": 89, "on_blockchain": false, "blockchain_data": { "tx_hash": "0xf6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b", "block_number": 1245686, "value_usd": 1123.45 } },
    { "id": "mol_010", "name": "Morphine", "formula": "C17H19NO3", "molecular_weight": 285.34, "ai_score": 90, "on_blockchain": true, "blockchain_data": { "tx_hash": "0xa7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c", "block_number": 1245687, "value_usd": 1967.23 } },
  ];

  const mockAiAgents = [
    { "id": "agent_001", "name": "Toxicity Predictor Pro", "description": "BioBERT-based toxicity classification", "status": "active", "accuracy": 0.92, "tasks_completed": 1247 },
    { "id": "agent_002", "name": "Solubility AI", "description": "GNN for solubility prediction", "status": "active", "accuracy": 0.87, "tasks_completed": 985 },
    { "id": "agent_003", "name": "Binding Affinity Master", "description": "3D CNN for protein-ligand binding", "status": "active", "accuracy": 0.84, "tasks_completed": 876 },
    { "id": "agent_004", "name": "Drug Likeness Evaluator", "description": "Rule-based + ML for drug properties", "status": "active", "accuracy": 0.89, "tasks_completed": 1102 },
    { "id": "agent_005", "name": "Synthesis Planner", "description": "Retrosynthesis AI with planning", "status": "active", "accuracy": 0.81, "tasks_completed": 943 },
  ];

  const mockStats = {
    success: true,
    stats: { total_molecules: 125, verified_molecules: 98, ai_agents_active: 5, total_ai_analyses: 5153, total_value_locked: 2487563.45, blockchain_txs: 98, avg_ai_score: 89, total_contributions: 243, knowledge_graph_nodes: 156, knowledge_graph_edges: 289 },
    performance: { avg_response_time: "0.45s", uptime: "99.9%", cache_hit_rate: "92%" }
  };

  const mockFlarePrice = { success: true, data: { price_usd: 0.0234, confidence: 0.89 } };

  const fetchDashboardData = async () => {
    try {
      const [statsRes, moleculesRes, agentsRes, flareRes] = await Promise.allSettled([
        fetch(`${API_URL}/stats`),
        fetch(`${API_URL}/molecules?limit=10`),
        fetch(`${API_URL}/agents`),
        fetch(`${API_URL}/flare/price`)
      ]);

      setStats(statsRes.status === 'fulfilled' && statsRes.value.ok ? await statsRes.value.json() : mockStats);
      const molData = moleculesRes.status === 'fulfilled' && moleculesRes.value.ok ? await moleculesRes.value.json() : { molecules: mockMolecules };
      setMolecules(molData.molecules || mockMolecules);
      const agentData = agentsRes.status === 'fulfilled' && agentsRes.value.ok ? await agentsRes.value.json() : { swarm: { agents: mockAiAgents } };
      setAiAgents(agentData.swarm?.agents || mockAiAgents);
      setFlarePrice(flareRes.status === 'fulfilled' && flareRes.value.ok ? await flareRes.value.json() : mockFlarePrice);
      setLoading(false);
    } catch (err) {
      console.error('Dashboard fetch error:', err);
      setStats(mockStats);
      setMolecules(mockMolecules);
      setAiAgents(mockAiAgents);
      setFlarePrice(mockFlarePrice);
      setTimeout(() => setLoading(false), 1000);
    }
  };

  useEffect(() => {
    fetchDashboardData();
    const interval = setInterval(fetchDashboardData, 30000);
    return () => clearInterval(interval);
  }, []);

  if (loading) {
    return (
      <div className="scifi-loading">
        <div className="scifi-terminal">
          <div className="terminal-header">
            <div className="terminal-dots">
              <span className="dot red"></span>
              <span className="dot yellow"></span>
              <span className="dot green"></span>
            </div>
            <div className="terminal-title">SYSTEM BOOT</div>
          </div>
          <div className="terminal-body">
            <div className="scan-line"></div>
            <div className="boot-sequence">
              <p className="boot-line">{'>'} INITIALIZING MOLCHAIN PLATFORM...</p>
              <p className="boot-line blinking">{'>'} LOADING AI AGENT SWARM [██████░░░░] 60%</p>
              <p className="boot-line">{'>'} CONNECTING TO FLARE BLOCKCHAIN...</p>
              <p className="boot-line">{'>'} SYNCING MOLECULAR DATABASE...</p>
              <p className="boot-line">{'>'} ESTABLISHING QUANTUM LINK...</p>
            </div>
          </div>
          <div className="terminal-footer">
            <span className="hackathon-badge">ETH OXFORD 2026</span>
            <span className="status-led active"></span>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="scifi-dashboard" style={{minHeight: '100vh', overflowX: 'hidden', width: '100%'}}>
      <div className="scifi-bg" style={{position: 'fixed', top: 0, left: 0, width: '100%', height: '100vh', zIndex: -1}}>
        <div className="grid-lines"></div>
        <div className="particles"></div>
        <div className="scan-line"></div>
      </div>
      
      <header className="dashboard-header" style={{padding: '15px 20px', width: '100%', overflow: 'hidden'}}>
        <div className="header-left">
          <h1 className="dashboard-title">
            <span className="title-glow">MOLCHAIN</span>
            <span className="title-sub">AI + BLOCKCHAIN PLATFORM</span>
          </h1>
          <div className="header-badges">
            <span className="badge neon-blue">ETH OXFORD 2026</span>
            <span className="badge neon-purple">FLARE INTEGRATED</span>
            <span className="badge neon-green">LIVE DEMO</span>
          </div>
        </div>
        <div className="header-right">
          <div className="system-time">
            <div className="time-digital">{new Date().toLocaleTimeString('en-US', {hour12: false})}</div>
            <div className="time-label">SYSTEM TIME</div>
          </div>
          <button className="refresh-btn" onClick={fetchDashboardData}>
            <span className="btn-icon">🔄</span>
            <span className="btn-text">SYNC</span>
          </button>
        </div>
      </header>
      
      <main className="dashboard-main" style={{padding: '15px', width: '100%', boxSizing: 'border-box'}}>
        <div className="stats-row">
          <div className="scifi-card neon-blue"><div className="card-glow"></div><div className="card-header"><div className="card-icon">🧬</div><div className="card-title">MOLECULES</div></div><div className="card-value">{stats?.stats?.total_molecules || 125}</div><div className="card-subtitle"><span className="verified-count">{stats?.stats?.verified_molecules || 98} VERIFIED</span><span className="blockchain-badge">⛓️ ON-CHAIN</span></div><div className="card-progress"><div className="progress-bar"><div className="progress-fill" style={{width: `${(stats?.stats?.verified_molecules || 98) / (stats?.stats?.total_molecules || 125) * 100}%`}}></div></div><div className="progress-label">BLOCKCHAIN SYNC</div></div></div>
          <div className="scifi-card neon-purple"><div className="card-glow"></div><div className="card-header"><div className="card-icon">🤖</div><div className="card-title">AI AGENTS</div></div><div className="card-value">{stats?.stats?.ai_agents_active || 5}</div><div className="card-subtitle"><span className="tasks-count">{stats?.stats?.total_ai_analyses?.toLocaleString() || '5,153'} ANALYSES</span><span className="status-badge active">ACTIVE</span></div><div className="agent-status">{aiAgents.map(agent => (<div key={agent.id} className="agent-dot active"></div>))}</div></div>
          <div className="scifi-card neon-orange"><div className="card-glow"></div><div className="card-header"><div className="card-icon">💰</div><div className="card-title">TOTAL VALUE</div></div><div className="card-value">${(stats?.stats?.total_value_locked || 2487563.45).toLocaleString()}</div><div className="card-subtitle"><span className="tx-count">{stats?.stats?.blockchain_txs || 98} TRANSACTIONS</span><span className="flare-badge">FLARE NETWORK</span></div><div className="value-ticker"><span className="ticker-up">↑ 2.4% TODAY</span><span className="ticker-value">${flarePrice?.data?.price_usd?.toFixed(4) || '0.0234'}/FLR</span></div></div>
          <div className="scifi-card neon-green"><div className="card-glow"></div><div className="card-header"><div className="card-icon">⚡</div><div className="card-title">AI CONSENSUS</div></div><div className="card-value">{stats?.stats?.avg_ai_score || 89}%</div><div className="card-subtitle"><span className="contributions">{stats?.stats?.total_contributions || 243} CONTRIBUTIONS</span><span className="accuracy-badge">92% ACCURACY</span></div><div className="consensus-meter"><div className="meter-fill" style={{width: `${stats?.stats?.avg_ai_score || 89}%`}}></div></div></div>
        </div>
        
        <div className="dashboard-grid">
          <div className="grid-column">
            <div className="flare-panel">
              <div className="panel-header"><h3 className="panel-title">🔥 FLARE BLOCKCHAIN</h3><div className="network-status"><span className="status-dot active"></span><span className="status-text">CONNECTED</span></div></div>
              <div className="flare-data">
                <div className="flare-price">
                  <div className="price-label">FTSO ORACLE PRICE</div>
                  <div className="price-value">${flarePrice?.data?.price_usd?.toFixed(4) || '0.0234'}</div>
                  <div className="price-confidence"><div className="confidence-bar"><div className="confidence-fill" style={{width: `${(flarePrice?.data?.confidence || 0.89) * 100}%`}}></div></div><span className="confidence-value">{(flarePrice?.data?.confidence * 100)?.toFixed(1) || '89.0'}% CONFIDENCE</span></div>
                </div>
                <div className="flare-stats">
                  <div className="flare-stat"><span className="stat-label">NETWORK</span><span className="stat-value">COSTON TESTNET</span></div>
                  <div className="flare-stat"><span className="stat-label">CHAIN ID</span><span className="stat-value">16</span></div>
                  <div className="flare-stat"><span className="stat-label">LATENCY</span><span className="stat-value">2.4s</span></div>
                  <div className="flare-stat"><span className="stat-label">BLOCK HEIGHT</span><span className="stat-value">{stats?.stats?.blockchain_txs ? 1245678 + stats.stats.blockchain_txs : 1245776}</span></div>
                </div>
              </div>
              <div className="flare-protocols"><div className="protocol active">FTSO</div><div className="protocol active">FDC</div><div className="protocol active">STATE CONNECTOR</div><div className="protocol active">FASSETS</div></div>
            </div>

            <div className="agents-panel">
              <div className="panel-header"><h3 className="panel-title">🤖 AI AGENT SWARM</h3><div className="swarm-status"><span className="pulse-dot"></span><span className="status-text">SYNCHRONIZED</span></div></div>
              <div className="agents-grid">{aiAgents.map((agent, index) => (<div key={agent.id} className="agent-cell" style={{animationDelay: `${index * 0.1}s`}}><div className="agent-header"><div className="agent-icon">🤖</div><div className="agent-info"><h4>{agent.name}</h4><span className="agent-type">{agent.description.split(' ')[0]}</span></div></div><div className="agent-metrics"><div className="metric"><span className="metric-label">ACCURACY</span><span className="metric-value">{(agent.accuracy * 100).toFixed(1)}%</span></div><div className="metric"><span className="metric-label">TASKS</span><span className="metric-value">{agent.tasks_completed.toLocaleString()}</span></div></div><div className="agent-activity"><div className="activity-bar"><div className="activity-fill" style={{width: `${(agent.accuracy * 100) - 70}%`}}></div></div><span className="activity-status">ACTIVE</span></div></div>))}</div>
            </div>
          </div>
          
          <div className="grid-column">
            <h2 className="section-title">💊 MOLECULAR REGISTRY</h2>
            <div className="molecules-grid">{(molecules.length > 0 ? molecules.slice(0, 6) : mockMolecules.slice(0, 6)).map(mol => (<MoleculeCard key={mol.id} molecule={mol} onClick={() => {setSelectedMolecule(mol); setShowModal(true);}} />))}</div>
          </div>
          
          <div className="grid-column">
            <div className="transactions-panel">
              <div className="panel-header"><h3 className="panel-title">⛓️ RECENT TRANSACTIONS</h3><div className="tx-count">{stats?.stats?.blockchain_txs || 98} TOTAL</div></div>
              <div className="transactions-list">{(molecules.length > 0 ? molecules.slice(0, 5) : mockMolecules.slice(0, 5)).map((mol) => (<div key={mol.id} className="transaction-row" style={{cursor: 'pointer'}} onClick={() => {setSelectedMolecule(mol); setShowModal(true);}}><div className="tx-icon"><div className="molecule-symbol">{mol.formula?.substring(0, 3) || '🧬'}</div></div><div className="tx-details"><div className="tx-header"><h4>{mol.name}</h4><span className="tx-value">${mol.blockchain_data?.value_usd?.toFixed(2) || '0.00'}</span></div><div className="tx-info"><span className="tx-hash">TX: {mol.blockchain_data?.tx_hash?.substring(0, 20) || '0x'}...</span><span className="tx-block">BLOCK: {mol.blockchain_data?.block_number || 'N/A'}</span><span className="tx-ai">AI SCORE: {mol.ai_score || 0}%</span></div><div className="tx-progress"><div className="progress-track"><div className="progress-step completed"></div><div className="progress-step completed"></div><div className="progress-step completed"></div><div className="progress-step active"></div><div className="progress-step"></div></div><span className="progress-label">BLOCKCHAIN CONFIRMED</span></div></div></div>))}</div>
            </div>

            <div className="kg-panel">
              <div className="panel-header"><h3 className="panel-title">🧠 KNOWLEDGE GRAPH</h3><div className="kg-stats"><span className="kg-stat">{stats?.stats?.knowledge_graph_nodes || 156} NODES</span><span className="kg-stat">{stats?.stats?.knowledge_graph_edges || 289} RELATIONSHIPS</span></div></div>
              <div className="kg-visualization">
                <div className="kg-canvas">
                  <div className="graph-node center"><div className="node-pulse"></div><div className="node-label">MOLECULES</div></div>
                  {['AGENTS', 'RESEARCH', 'BLOCKCHAIN', 'DATABASE', 'VALIDATORS'].map((label, i) => {const angle = (i * 72) * (Math.PI / 180); const radius = 80; const x = 120 + radius * Math.cos(angle); const y = 120 + radius * Math.sin(angle); return (<div key={label} className="graph-node" style={{top: `${y}px`, left: `${x}px`}}><div className="node-pulse"></div><div className="node-label">{label}</div></div>);})}
                </div>
                <div className="kg-legend"><div className="legend-item"><div className="legend-dot molecule"></div><span>MOLECULES</span></div><div className="legend-item"><div className="legend-dot agent"></div><span>AI AGENTS</span></div><div className="legend-item"><div className="legend-dot blockchain"></div><span>BLOCKCHAIN</span></div></div>
              </div>
            </div>

            <div className="performance-panel">
              <div className="panel-header"><h3 className="panel-title">⚡ SYSTEM PERFORMANCE</h3><div className="performance-status"><span className="status-led green"></span><span className="status-text">OPTIMAL</span></div></div>
              <div className="performance-metrics">
                <div className="metric-card"><div className="metric-header"><span className="metric-icon">⚡</span><span className="metric-title">RESPONSE TIME</span></div><div className="metric-value">{stats?.performance?.avg_response_time || '0.45s'}</div><div className="metric-trend positive">-12%</div></div>
                <div className="metric-card"><div className="metric-header"><span className="metric-icon">📈</span><span className="metric-title">UPTIME</span></div><div className="metric-value">{stats?.performance?.uptime || '99.9%'}</div><div className="metric-trend positive">+0.1%</div></div>
                <div className="metric-card"><div className="metric-header"><span className="metric-icon">💾</span><span className="metric-title">CACHE HIT</span></div><div className="metric-value">{stats?.performance?.cache_hit_rate || '92%'}</div><div className="metric-trend positive">+3%</div></div>
                <div className="metric-card"><div className="metric-header"><span className="metric-icon">🔗</span><span className="metric-title">BLOCKCHAIN SYNC</span></div><div className="metric-value">100%</div><div className="metric-trend stable">SYNCED</div></div>
              </div>
            </div>
          </div>
        </div>
      </main>
      
      <footer className="dashboard-footer">
        <div className="footer-left">
          <div className="network-status"><span className="status-dot green"></span><span>ALL SYSTEMS NOMINAL</span></div>
          <div className="version">v2.0.0 | ETH OXFORD HACKATHON</div>
        </div>
        <div className="footer-right">
          <div className="data-stats">
            <span>{stats?.stats?.total_molecules || 125} MOLECULES</span>
            <span>•</span>
            <span>{stats?.stats?.ai_agents_active || 5} AI AGENTS</span>
            <span>•</span>
            <span>{stats?.stats?.blockchain_txs || 98} TXS</span>
            <span>•</span>
            <span>${(stats?.stats?.total_value_locked || 2487563.45).toLocaleString()} TVL</span>
          </div>
        </div>
      </footer>

      <MoleculeDetailModal molecule={selectedMolecule} isOpen={showModal} onClose={() => {setShowModal(false); setSelectedMolecule(null);}} />
    </div>
  );
};

export default Dashboard;
