from flask import Flask, jsonify, request
from flask_cors import CORS
import random
import time
from datetime import datetime

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": ["http://localhost:5173", "http://localhost:5173", "http://127.0.0.1:5173"]}})  # Enable CORS for React frontend

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response


# ===== IN-MEMORY DATABASES =====
molecules_db = []
contributions_db = []

# ===== INITIAL SAMPLE DATA =====
def init_sample_data():
    molecules_db.extend([
        {
            "id": "aspirin_001",
            "name": "Aspirin",
            "formula": "C9H8O4",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "weight": 180.16,
            "description": "Pain reliever, anti-inflammatory, antiplatelet",
            "ai_score": 87,
            "on_blockchain": True,
            "tx_hash": "0x7a3f...b89c",
            "block_number": 1245678,
            "properties": {
                "toxicity": "Low",
                "solubility": "Moderate",
                "melting_point": "135°C",
                "bioavailability": "High"
            },
            "contributions_count": 3,
            "created_at": "2024-01-15T10:30:00Z",
            "researcher": "Dr. Jane Smith"
        },
        {
            "id": "caffeine_001",
            "name": "Caffeine",
            "formula": "C8H10N4O2",
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "weight": 194.19,
            "description": "Central nervous system stimulant",
            "ai_score": 92,
            "on_blockchain": True,
            "tx_hash": "0x8b4e...c67d",
            "block_number": 1245679,
            "properties": {
                "toxicity": "Moderate",
                "solubility": "Moderate",
                "melting_point": "238°C",
                "bioavailability": "High"
            },
            "contributions_count": 5,
            "created_at": "2024-01-16T14:20:00Z",
            "researcher": "Dr. Robert Chen"
        },
        {
            "id": "ibuprofen_001",
            "name": "Ibuprofen",
            "formula": "C13H18O2",
            "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(O)=O",
            "weight": 206.28,
            "description": "Nonsteroidal anti-inflammatory drug",
            "ai_score": 78,
            "on_blockchain": True,
            "tx_hash": "0x9c5f...d78e",
            "block_number": 1245680,
            "properties": {
                "toxicity": "Low",
                "solubility": "Low",
                "melting_point": "76°C",
                "bioavailability": "Moderate"
            },
            "contributions_count": 2,
            "created_at": "2024-01-17T09:15:00Z",
            "researcher": "Dr. Maria Garcia"
        },
        {
            "id": "vitamin_c_001",
            "name": "Vitamin C",
            "formula": "C6H8O6",
            "smiles": "C(C(C1C(=C(C(=O)O1)O)O)O)O",
            "weight": 176.12,
            "description": "Essential nutrient, antioxidant",
            "ai_score": 95,
            "on_blockchain": False,
            "tx_hash": None,
            "block_number": None,
            "properties": {
                "toxicity": "Very Low",
                "solubility": "High",
                "melting_point": "190°C",
                "bioavailability": "High"
            },
            "contributions_count": 0,
            "created_at": "2024-01-18T11:45:00Z",
            "researcher": "Dr. Alan Turing"
        }
    ])
    
    contributions_db.extend([
        {
            "id": "cont_001",
            "molecule_id": "aspirin_001",
            "title": "Quantum Mechanical Analysis",
            "type": "analysis",
            "researcher": "Dr. Alan Turing",
            "description": "DFT calculations reveal new binding sites",
            "methodology": "Density Functional Theory, B3LYP/6-31G*",
            "results": "Identified 3 novel binding conformations",
            "data_url": "https://example.com/data1.zip",
            "ai_validation": 94,
            "timestamp": "2024-01-20T11:30:00Z"
        },
        {
            "id": "cont_002",
            "molecule_id": "caffeine_001",
            "title": "Solubility Enhancement Study",
            "type": "experiment",
            "researcher": "Dr. Marie Curie",
            "description": "Experimental study on solubility in various solvents",
            "methodology": "UV-Vis spectroscopy, HPLC",
            "results": "Found optimal solubility in DMSO",
            "data_url": "https://example.com/data2.zip",
            "ai_validation": 88,
            "timestamp": "2024-01-19T14:20:00Z"
        }
    ])

init_sample_data()

# ===== HEALTH CHECK =====
@app.route('/api/health', methods=['GET'])
def health():
    return jsonify({
        "status": "healthy",
        "service": "MolChain AI Backend",
        "version": "1.0.0",
        "timestamp": datetime.now().isoformat() + "Z",
        "stats": {
            "molecules": len(molecules_db),
            "contributions": len(contributions_db),
            "users": 0
        }
    })

# ===== GET ALL MOLECULES =====
@app.route('/api/molecules', methods=['GET'])
def get_molecules():
    sort_by = request.args.get('sort', 'ai_score')
    filter_type = request.args.get('filter', 'all')
    
    filtered = molecules_db.copy()
    
    if filter_type == 'verified':
        filtered = [m for m in filtered if m['on_blockchain']]
    elif filter_type == 'pending':
        filtered = [m for m in filtered if not m['on_blockchain']]
    
    if sort_by == 'ai_score':
        filtered.sort(key=lambda x: x['ai_score'], reverse=True)
    elif sort_by == 'name':
        filtered.sort(key=lambda x: x['name'])
    elif sort_by == 'recent':
        filtered.sort(key=lambda x: x.get('created_at', ''), reverse=True)
    
    return jsonify({
        "success": True,
        "count": len(filtered),
        "molecules": filtered
    })

# ===== GET SINGLE MOLECULE =====
@app.route('/api/molecules/<molecule_id>', methods=['GET'])
def get_molecule(molecule_id):
    molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
    if molecule:
        # Get contributions for this molecule
        molecule_contributions = [c for c in contributions_db if c["molecule_id"] == molecule_id]
        
        response = {
            "success": True,
            "molecule": molecule,
            "contributions": molecule_contributions,
            "ai_agents": [
                {"name": "Toxicity AI", "status": "active", "last_analysis": "2024-01-20T10:30:00Z"},
                {"name": "Solubility Predictor", "status": "active", "last_analysis": "2024-01-20T10:31:00Z"},
                {"name": "Binding Affinity AI", "status": "active", "last_analysis": "2024-01-20T10:32:00Z"},
                {"name": "Stability Checker", "status": "active", "last_analysis": "2024-01-20T10:33:00Z"}
            ]
        }
        return jsonify(response)
    return jsonify({"success": False, "error": "Molecule not found"}), 404

# ===== CREATE NEW MOLECULE =====
@app.route('/api/molecules', methods=['POST'])
def create_molecule():
    data = request.json
    if not data.get('name') or not data.get('formula'):
        return jsonify({"success": False, "error": "Missing required fields"}), 400
    
    molecule_id = f"{data['name'].lower().replace(' ', '_')}_{random.randint(1000, 9999)}"
    
    new_molecule = {
        "id": molecule_id,
        "name": data['name'],
        "formula": data['formula'],
        "smiles": data.get('smiles', ''),
        "weight": data.get('weight', 0),
        "description": data.get('description', ''),
        "ai_score": 0,
        "on_blockchain": False,
        "tx_hash": None,
        "block_number": None,
        "properties": {
            "toxicity": "Pending",
            "solubility": "Pending",
            "melting_point": "Pending",
            "bioavailability": "Pending"
        },
        "contributions_count": 0,
        "created_at": datetime.now().isoformat() + "Z",
        "researcher": data.get('researcher', 'Anonymous')
    }
    
    molecules_db.append(new_molecule)
    
    return jsonify({
        "success": True,
        "message": "Molecule created successfully",
        "molecule": new_molecule,
        "next_step": "Run AI analysis to validate properties"
    }), 201

# ===== AI ANALYSIS ENDPOINT (CRITICAL FOR REACT) =====
@app.route('/api/analyze', methods=['POST'])
def analyze_molecule():
    data = request.json
    molecule_id = data.get('molecule_id')
    
    if not molecule_id:
        return jsonify({"success": False, "error": "Missing molecule_id"}), 400
    
    molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
    if not molecule:
        return jsonify({"success": False, "error": "Molecule not found"}), 404
    
    # Simulate AI agent swarm analysis
    print(f"🤖 Starting AI analysis for {molecule['name']}...")
    
    agents = [
        {"name": "Toxicity AI", "task": "toxicity", "delay": 0.5},
        {"name": "Solubility Predictor", "task": "solubility", "delay": 0.7},
        {"name": "Binding Affinity AI", "task": "binding", "delay": 0.9},
        {"name": "Stability Checker", "task": "stability", "delay": 1.1},
        {"name": "Bioavailability AI", "task": "bioavailability", "delay": 1.3}
    ]
    
    results = []
    total_weighted_score = 0
    total_weight = 0
    
    for agent in agents:
        time.sleep(agent['delay'])
        score = random.randint(70, 98)
        confidence = round(random.uniform(0.8, 0.98), 2)
        
        agent_result = {
            "agent": agent["name"],
            "task": agent["task"],
            "score": score,
            "confidence": confidence,
            "status": "✅ Valid" if score > 75 else "⚠️ Warning",
            "details": f"Analyzed {agent['task']} using quantum mechanical simulations"
        }
        
        results.append(agent_result)
        total_weighted_score += score * 0.2  # Equal weight for all agents
        total_weight += 0.2
    
    consensus_score = round(total_weighted_score / total_weight)
    
    # Update molecule with AI results
    molecule['ai_score'] = consensus_score
    molecule['on_blockchain'] = True
    molecule['tx_hash'] = f"0x{random.getrandbits(256):064x}"[:20] + "..." + f"{random.getrandbits(256):064x}"[-4:]
    molecule['block_number'] = random.randint(1200000, 1300000)
    molecule['properties'] = {
        "toxicity": "Low" if results[0]['score'] > 80 else "Moderate",
        "solubility": "High" if results[1]['score'] > 85 else "Moderate",
        "melting_point": f"{random.randint(50, 300)}°C",
        "bioavailability": "High" if results[4]['score'] > 85 else "Moderate"
    }
    
    print(f"✅ AI analysis complete. Consensus score: {consensus_score}%")
    
    return jsonify({
        "success": True,
        "message": f"AI analysis complete for {molecule['name']}",
        "consensus_score": consensus_score,
        "agents": results,
        "blockchain": {
            "tx_hash": molecule['tx_hash'],
            "block_number": molecule['block_number'],
            "status": "confirmed",
            "explorer_url": f"https://flare-explorer.example.com/tx/{molecule['tx_hash']}"
        },
        "molecule": molecule
    })

# ===== STATS ENDPOINT (CRITICAL FOR REACT) =====
@app.route('/api/stats', methods=['GET'])
def get_stats():
    total_ai_score = sum(m['ai_score'] for m in molecules_db if m['ai_score'] > 0)
    avg_ai_score = round(total_ai_score / len(molecules_db)) if molecules_db else 0
    
    return jsonify({
        "success": True,
        "stats": {
            "total_molecules": len(molecules_db),
            "verified_molecules": len([m for m in molecules_db if m['on_blockchain']]),
            "total_contributions": len(contributions_db),
            "avg_ai_score": avg_ai_score,
            "ai_agents_active": 5,
            "total_analyses": len(molecules_db) * 5,
            "blockchain_txs": len([m for m in molecules_db if m['tx_hash']])
        },
        "timestamp": datetime.now().isoformat() + "Z"
    })

# ===== AI AGENTS ENDPOINT (CRITICAL FOR REACT) =====
@app.route('/api/agents', methods=['GET'])
def get_agents():
    return jsonify({
        "success": True,
        "swarm": {
            "name": "MolChain AI Agent Swarm v1.0",
            "status": "active",
            "agents": [
                {"name": "Toxicity AI", "status": "active", "accuracy": 0.92, "tasks": 1247, "last_active": "2024-01-20T10:30:00Z"},
                {"name": "Solubility Predictor", "status": "active", "accuracy": 0.87, "tasks": 985, "last_active": "2024-01-20T10:31:00Z"},
                {"name": "Binding Affinity AI", "status": "active", "accuracy": 0.84, "tasks": 876, "last_active": "2024-01-20T10:32:00Z"},
                {"name": "Stability Checker", "status": "active", "accuracy": 0.91, "tasks": 1102, "last_active": "2024-01-20T10:33:00Z"},
                {"name": "Bioavailability AI", "status": "active", "accuracy": 0.88, "tasks": 943, "last_active": "2024-01-20T10:34:00Z"}
            ],
            "total_tasks": 5153,
            "consensus_accuracy": 0.89
        }
    })

# ===== FLARE PRICE ENDPOINT (CRITICAL FOR REACT) =====
@app.route('/api/flare/price', methods=['GET'])
def get_flare_price():
    """Simulate Flare FTSO price feed"""
    return jsonify({
        "success": True,
        "oracle": "Flare FTSO",
        "price_usd": round(random.uniform(0.02, 0.05), 4),
        "confidence": round(random.uniform(0.85, 0.99), 2),
        "timestamp": datetime.now().isoformat() + "Z"
    })

# ===== FLARE MOLECULE VALUE ENDPOINT =====
@app.route('/api/flare/<molecule_id>', methods=['GET'])
def get_molecule_value(molecule_id):
    """Get molecular valuation from Flare oracle"""
    molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
    if not molecule:
        return jsonify({"success": False, "error": "Molecule not found"}), 404
    
    base_value = molecule['ai_score'] * random.uniform(100, 500)
    
    return jsonify({
        "success": True,
        "molecule": molecule_id,
        "value_usd": round(base_value, 2),
        "valuation_model": "AI Consensus × Market Factor",
        "oracle": "Flare FTSO Enhanced",
        "confidence": round(random.uniform(0.7, 0.95), 2),
        "timestamp": datetime.now().isoformat() + "Z"
    })

# ===== CONTRIBUTIONS ENDPOINTS =====
@app.route('/api/contributions', methods=['GET'])
def get_contributions():
    molecule_id = request.args.get('molecule_id')
    
    if molecule_id:
        filtered = [c for c in contributions_db if c["molecule_id"] == molecule_id]
    else:
        filtered = contributions_db
    
    return jsonify({
        "success": True,
        "count": len(filtered),
        "contributions": filtered
    })

@app.route('/api/contributions', methods=['POST'])
def create_contribution():
    data = request.json
    
    required_fields = ['molecule_id', 'title', 'type', 'researcher']
    for field in required_fields:
        if not data.get(field):
            return jsonify({"success": False, "error": f"Missing field: {field}"}), 400
    
    contribution_id = f"cont_{len(contributions_db) + 1:03d}"
    
    new_contribution = {
        "id": contribution_id,
        "molecule_id": data['molecule_id'],
        "title": data['title'],
        "type": data['type'],
        "researcher": data['researcher'],
        "description": data.get('description', ''),
        "methodology": data.get('methodology', ''),
        "results": data.get('results', ''),
        "data_url": data.get('data_url', ''),
        "ai_validation": random.randint(70, 98),
        "timestamp": datetime.now().isoformat() + "Z"
    }
    
    contributions_db.append(new_contribution)
    
    # Update molecule contribution count
    for molecule in molecules_db:
        if molecule['id'] == data['molecule_id']:
            molecule['contributions_count'] = molecule.get('contributions_count', 0) + 1
            break
    
    return jsonify({
        "success": True,
        "message": "Contribution submitted successfully",
        "contribution": new_contribution
    }), 201

# ===== ROOT ENDPOINT =====
@app.route('/')
def home():
    return """
    <html>
    <head>
        <title>MolChain Backend</title>
        <style>
            body { 
                background: #0f172a; 
                color: white; 
                font-family: sans-serif; 
                padding: 40px;
                max-width: 800px;
                margin: 0 auto;
            }
            h1 { color: #60a5fa; }
            .endpoint { 
                background: #1e293b; 
                padding: 15px; 
                margin: 10px 0; 
                border-radius: 8px;
                border-left: 4px solid #3b82f6;
            }
            .method { 
                display: inline-block; 
                padding: 4px 8px; 
                background: #3b82f6; 
                color: white; 
                border-radius: 4px;
                font-weight: bold;
                margin-right: 10px;
            }
            .url { color: #cbd5e1; font-family: monospace; }
        </style>
    </head>
    <body>
        <h1>🧬 MolChain Python Backend</h1>
        <p>Backend server for the MolChain AI + Blockchain platform</p>
        
        <h2>📡 Available Endpoints:</h2>
        
        <div class="endpoint">
            <span class="method">GET</span>
            <span class="url">/api/health</span>
            <p>Health check and service status</p>
        </div>
        
        <div class="endpoint">
            <span class="method">GET</span>
            <span class="url">/api/molecules</span>
            <p>Get all molecules (use ?sort=ai_score&filter=verified)</p>
        </div>
        
        <div class="endpoint">
            <span class="method">POST</span>
            <span class="url">/api/analyze</span>
            <p>Run AI analysis on a molecule (send {"molecule_id": "aspirin_001"})</p>
        </div>
        
        <div class="endpoint">
            <span class="method">GET</span>
            <span class="url">/api/stats</span>
            <p>Get platform statistics</p>
        </div>
        
        <div class="endpoint">
            <span class="method">GET</span>
            <span class="url">/api/agents</span>
            <p>Get AI agent swarm status</p>
        </div>
        
        <div class="endpoint">
            <span class="method">GET</span>
            <span class="url">/api/flare/price</span>
            <p>Get simulated Flare FTSO price feed</p>
        </div>
        
        <p>🚀 Frontend should connect to: <code>http://localhost:5000/api</code></p>
    </body>
    </html>
    """


# ===== START SERVER =====
if __name__ == '__main__':
    print("\n" + "="*60)
    print("🧬 MOLCHAIN AI BACKEND STARTING...")
    print("="*60)
    print("\n📡 Essential Endpoints:")
    print("  Health:      http://localhost:5173/api/health")
    print("  Molecules:   http://localhost:5173/api/molecules")
    print("  Analyze:     POST http://localhost:5173/api/analyze")
    print("  Stats:       http://localhost:5173/api/stats")
    print("  Agents:      http://localhost:5173/api/agents")
    print("  Flare Price: http://localhost:5173/api/flare/price")
    print("\n" + "="*60)
    print("✅ Server ready! React frontend should connect automatically")
    print("⏰ Press Ctrl+C to stop the server")
    print("="*60 + "\n")
    
    app.run(debug=True, port=5173, use_reloader=False)