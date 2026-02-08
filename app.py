from flask import Flask, jsonify, request, send_file
from flask_cors import CORS
import random
import time
import json
import os
from datetime import datetime, timedelta
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import wraps
import hashlib
import base64
from io import BytesIO
from molecules_database import database


# AI/ML imports
import torch
from transformers import pipeline
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import threading

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": ["http://localhost:5173", "http://localhost:3000", "http://127.0.0.1:5173", "*"]}})

# ===== ENHANCED CONFIGURATION =====
class Config:
    AI_AGENTS = 5
    FLARE_TESTNET = "https://coston-api.flare.network/ext/C/rpc"
    MAX_WORKERS = 10
    SESSION_TIMEOUT = 300
    MAX_MOLECULES = 1000
    VERSION = "2.0.0-HACKATHON"

app.config.from_object(Config)

# ===== RATE LIMITING DECORATOR =====
def rate_limit(max_per_minute=60):
    def decorator(f):
        calls = []
        @wraps(f)
        def wrapper(*args, **kwargs):
            now = time.time()
            calls.append(now)
            calls[:] = [call for call in calls if now - call < 60]
            if len(calls) > max_per_minute:
                return jsonify({"error": "Rate limit exceeded"}), 429
            return f(*args, **kwargs)
        return wrapper
    return decorator

# ===== IN-MEMORY DATABASES WITH ENHANCEMENTS =====
molecules_db = []
contributions_db = []
ai_agents_db = []
knowledge_graph = []
blockchain_txs = []
web_scraping_jobs = []
researchers_db = []


# Lightweight `database` helper to provide simple accessors used by frontend
class Database:
    @property
    def molecules(self):
        return molecules_db

    @property
    def knowledge_graph(self):
        # If KnowledgeGraph object exists, expose nodes/edges, otherwise fallback
        try:
            return {
                'nodes': knowledge_graph.nodes,
                'edges': knowledge_graph.edges
            }
        except Exception:
            return {
                'nodes': knowledge_graph if isinstance(knowledge_graph, list) else [],
                'edges': []
            }

    def get_molecule_ecosystem(self, molecule_id):
        mol = next((m for m in molecules_db if m.get('id') == molecule_id), None)
        if not mol:
            return None

        contribs = [c for c in contributions_db if c.get('molecule_id') == molecule_id]
        agents = [a for a in ai_agents_db]
        kg_similar = []
        try:
            kg_similar = knowledge_graph.query('find_similar', molecule_id=molecule_id)
        except Exception:
            # best-effort: empty
            kg_similar = []

        return {
            'molecule': mol,
            'contributions': contribs,
            'agents': agents,
            'knowledge_graph_similar': kg_similar,
            'blockchain': mol.get('blockchain_data')
        }

    def search(self, query, category='all'):
        q = (query or '').strip().lower()
        results = []
        if category in ('all', 'molecule'):
            for m in molecules_db:
                if q in (m.get('name') or '').lower() or q in (m.get('smiles') or '').lower() or q in (m.get('formula') or '').lower():
                    results.append({'type': 'molecule', 'id': m.get('id'), 'name': m.get('name')})
        if category in ('all', 'researcher'):
            for r in researchers_db:
                if q in (r.get('name') or '').lower() or q in (r.get('email') or '').lower():
                    results.append({'type': 'researcher', 'id': r.get('id'), 'name': r.get('name')})
        return results

    def get_statistics(self):
        total_ai_score = sum(m.get('ai_score', 0) for m in molecules_db)
        avg_ai_score = round(total_ai_score / len(molecules_db)) if molecules_db else 0
        total_value = sum(m.get('blockchain_data', {}).get('value_usd', 0) for m in molecules_db if m.get('on_blockchain'))
        return {
            'total_molecules': len(molecules_db),
            'verified_molecules': len([m for m in molecules_db if m.get('on_blockchain')]),
            'total_contributions': len(contributions_db),
            'total_researchers': len(researchers_db),
            'avg_ai_score': avg_ai_score,
            'ai_agents_active': len([a for a in ai_agents_db if a.get('status') == 'active']),
            'total_value_locked': round(total_value, 2),
        }


database = Database()

# ===== AI AGENT SWARM IMPLEMENTATION =====
class AIAgentSwarm:
    def __init__(self):
        self.agents = {}
        self.executor = ThreadPoolExecutor(max_workers=app.config['AI_AGENTS'])
        self.initialize_agents()
        self.consensus_history = []
        
    def initialize_agents(self):
        """Initialize 5 specialized AI agents"""
        self.agents = {
            'toxicity': {
                'name': 'Toxicity Predictor',
                'model': self.load_toxicity_model(),
                'weight': 0.25,
                'accuracy': 0.92,
                'description': 'BioBERT-based toxicity classification',
                'tasks_completed': 0
            },
            'solubility': {
                'name': 'Solubility Predictor',
                'model': self.load_solubility_model(),
                'weight': 0.20,
                'accuracy': 0.87,
                'description': 'Graph neural network for solubility',
                'tasks_completed': 0
            },
            'binding': {
                'name': 'Binding Affinity AI',
                'model': self.load_binding_model(),
                'weight': 0.20,
                'accuracy': 0.84,
                'description': '3D CNN for protein-ligand interactions',
                'tasks_completed': 0
            },
            'druglikeness': {
                'name': 'Drug Likeness Evaluator',
                'model': self.load_druglikeness_model(),
                'weight': 0.20,
                'accuracy': 0.89,
                'description': 'Rule-based + ML for drug-like properties',
                'tasks_completed': 0
            },
            'synthesis': {
                'name': 'Synthesis Feasibility',
                'model': self.load_synthesis_model(),
                'weight': 0.15,
                'accuracy': 0.81,
                'description': 'Retrosynthesis planning AI',
                'tasks_completed': 0
            }
        }
        
        # Initialize agent database
        for agent_id, agent in self.agents.items():
            ai_agents_db.append({
                'id': agent_id,
                'name': agent['name'],
                'status': 'active',
                'accuracy': agent['accuracy'],
                'tasks_completed': 0,
                'last_active': datetime.now().isoformat() + "Z",
                'description': agent['description'],
                'weight': agent['weight']
            })
    
    def load_toxicity_model(self):
        """Load toxicity prediction model (FAKE/MOCKED)"""
        class ToxicityModel:
            def predict(self, smiles):
                # Return random score between 0.2 and 0.7
                return round(random.uniform(0.2, 0.7), 3)
        return ToxicityModel()
    
    def load_solubility_model(self):
        """Load solubility prediction model (FAKE/MOCKED)"""
        class SolubilityModel:
            def predict(self, smiles):
                # Return random score between 0.4 and 0.9
                return round(random.uniform(0.4, 0.9), 3)
        return SolubilityModel()
    
    def load_binding_model(self):
        """Load binding affinity model (FAKE/MOCKED)"""
        class BindingModel:
            def predict(self, smiles):
                # Return random score between 0.5 and 0.9
                return round(random.uniform(0.5, 0.9), 3)
        return BindingModel()
    
    def load_druglikeness_model(self):
        """Load drug-likeness model (FAKE/MOCKED)"""
        class DrugLikenessModel:
            def predict(self, smiles):
                # Return random score between 0.6 and 1.0
                return round(random.uniform(0.6, 1.0), 3)
        return DrugLikenessModel()
    
    def load_synthesis_model(self):
        """Load synthesis feasibility model (FAKE/MOCKED)"""
        class SynthesisModel:
            def predict(self, smiles):
                # Return random score between 0.5 and 0.95
                return round(random.uniform(0.5, 0.95), 3)
        return SynthesisModel()
    
    def analyze_molecule(self, smiles, molecule_id):
        """Run all agents in parallel and return consensus"""
        print(f"🤖 Starting AI swarm analysis for {smiles}")
        
        futures = {}
        for agent_id, agent in self.agents.items():
            future = self.executor.submit(
                self.run_agent_analysis,
                agent_id, agent, smiles, molecule_id
            )
            futures[future] = agent_id
        
        results = {}
        for future in as_completed(futures):
            agent_id = futures[future]
            try:
                result = future.result(timeout=10)
                results[agent_id] = result
                self.agents[agent_id]['tasks_completed'] += 1
            except Exception as e:
                results[agent_id] = {
                    'score': 0.5,
                    'confidence': 0.1,
                    'error': str(e)
                }
        
        # Calculate consensus
        consensus = self.calculate_consensus(results)
        
        # Store in consensus history
        self.consensus_history.append({
            'molecule_id': molecule_id,
            'smiles': smiles,
            'consensus': consensus,
            'timestamp': datetime.now().isoformat() + "Z",
            'agent_results': results
        })
        
        return {
            'agent_results': results,
            'consensus_score': consensus,
            'blockchain_ready': consensus > 0.7,
            'timestamp': datetime.now().isoformat() + "Z"
        }
    
    def run_agent_analysis(self, agent_id, agent, smiles, molecule_id):
        """Run individual agent analysis"""
        try:
            score = agent['model'].predict(smiles)
            confidence = agent['accuracy'] * random.uniform(0.8, 1.0)
            
            return {
                'agent_id': agent_id,
                'agent_name': agent['name'],
                'score': round(float(score), 3),
                'confidence': round(float(confidence), 3),
                'description': agent['description'],
                'timestamp': datetime.now().isoformat() + "Z",
                'status': 'success'
            }
        except Exception as e:
            return {
                'agent_id': agent_id,
                'agent_name': agent['name'],
                'score': 0.5,
                'confidence': 0.1,
                'error': str(e),
                'timestamp': datetime.now().isoformat() + "Z",
                'status': 'failed'
            }
    
    def calculate_consensus(self, results):
        """Calculate weighted consensus between agents"""
        total_score = 0
        total_weight = 0
        
        for agent_id, result in results.items():
            if 'error' not in result:
                weight = self.agents[agent_id]['weight']
                score = result['score']
                confidence = result.get('confidence', 0.5)
                
                total_score += score * weight * confidence
                total_weight += weight * confidence
        
        return round(total_score / total_weight if total_weight > 0 else 0.5, 3)

# Initialize AI Swarm
ai_swarm = AIAgentSwarm()

# ===== FLARE BLOCKCHAIN SIMULATION =====
class FlareBlockchainSimulator:
    def __init__(self):
        self.chain_id = 16  # Flare Coston
        self.gas_price = 25000000000  # 25 Gwei
        self.ftso_price = 0.023
        self.transactions = []
        self.molecule_registry = {}
        
    def register_molecule(self, molecule_data, ai_scores):
        """Simulate molecule registration on Flare blockchain"""
        tx_hash = f"0x{hashlib.sha256(str(molecule_data).encode()).hexdigest()[:64]}"
        block_number = random.randint(1500000, 1600000)
        
        # Calculate molecular value using FTSO
        base_value = ai_scores['consensus_score'] * 1000
        market_factor = random.uniform(0.8, 1.2)
        molecule_value = base_value * market_factor * self.ftso_price
        
        tx = {
            'tx_hash': tx_hash,
            'block_number': block_number,
            'molecule_id': molecule_data.get('id'),
            'status': 'confirmed',
            'timestamp': datetime.now().isoformat() + "Z",
            'value_usd': round(molecule_value, 2),
            'gas_used': random.randint(50000, 200000),
            'from': '0x' + ''.join(random.choices('0123456789abcdef', k=40)),
            'to': '0x' + ''.join(random.choices('0123456789abcdef', k=40)),
            'explorer_url': f"https://coston-explorer.flare.network/tx/{tx_hash}"
        }
        
        self.transactions.append(tx)
        self.molecule_registry[molecule_data.get('id')] = {
            'on_chain': True,
            'tx_hash': tx_hash,
            'block_number': block_number,
            'value_usd': tx['value_usd'],
            'registered_at': tx['timestamp']
        }
        
        return tx
    
    def get_ftso_price(self):
        """Simulate FTSO price feed"""
        # Add some random variation
        variation = random.uniform(-0.001, 0.001)
        self.ftso_price = max(0.01, self.ftso_price + variation)
        
        return {
            'price_usd': round(self.ftso_price, 4),
            'confidence': round(random.uniform(0.85, 0.99), 3),
            'timestamp': datetime.now().isoformat() + "Z",
            'oracle': 'FTSO v2',
            'decimals': 18
        }
    
    def get_molecule_value(self, molecule_id, ai_score):
        """Get molecular valuation from blockchain"""
        if molecule_id in self.molecule_registry:
            return self.molecule_registry[molecule_id]['value_usd']
        
        # Calculate if not registered
        base_value = ai_score * 1000
        market_factor = random.uniform(0.8, 1.2)
        return round(base_value * market_factor * self.ftso_price, 2)

flare_simulator = FlareBlockchainSimulator()

# ===== KNOWLEDGE GRAPH =====
class KnowledgeGraph:
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.initialize_sample_graph()
    
    def initialize_sample_graph(self):
        """Initialize sample knowledge graph for demo"""
        # Add molecules as nodes
        molecule_nodes = [
            {'id': 'aspirin', 'type': 'molecule', 'label': 'Aspirin', 'properties': {'formula': 'C9H8O4'}},
            {'id': 'caffeine', 'type': 'molecule', 'label': 'Caffeine', 'properties': {'formula': 'C8H10N4O2'}},
            {'id': 'ibuprofen', 'type': 'molecule', 'label': 'Ibuprofen', 'properties': {'formula': 'C13H18O2'}},
            {'id': 'paracetamol', 'type': 'molecule', 'label': 'Paracetamol', 'properties': {'formula': 'C8H9NO2'}},
            {'id': 'vitamin_c', 'type': 'molecule', 'label': 'Vitamin C', 'properties': {'formula': 'C6H8O6'}},
        ]
        
        # Add properties as nodes
        property_nodes = [
            {'id': 'pain_relief', 'type': 'property', 'label': 'Pain Relief'},
            {'id': 'anti_inflammatory', 'type': 'property', 'label': 'Anti-inflammatory'},
            {'id': 'stimulant', 'type': 'property', 'label': 'Stimulant'},
            {'id': 'antioxidant', 'type': 'property', 'label': 'Antioxidant'},
        ]
        
        # Add researchers as nodes
        researcher_nodes = [
            {'id': 'res_1', 'type': 'researcher', 'label': 'Dr. Jane Smith', 'properties': {'institution': 'Oxford'}},
            {'id': 'res_2', 'type': 'researcher', 'label': 'Dr. Robert Chen', 'properties': {'institution': 'Cambridge'}},
            {'id': 'res_3', 'type': 'researcher', 'label': 'Dr. Maria Garcia', 'properties': {'institution': 'ETH Zurich'}},
        ]
        
        self.nodes.extend(molecule_nodes + property_nodes + researcher_nodes)
        
        # Add relationships
        self.edges.extend([
            {'from': 'aspirin', 'to': 'pain_relief', 'relationship': 'HAS_PROPERTY', 'weight': 0.9},
            {'from': 'aspirin', 'to': 'anti_inflammatory', 'relationship': 'HAS_PROPERTY', 'weight': 0.8},
            {'from': 'caffeine', 'to': 'stimulant', 'relationship': 'HAS_PROPERTY', 'weight': 0.95},
            {'from': 'vitamin_c', 'to': 'antioxidant', 'relationship': 'HAS_PROPERTY', 'weight': 0.98},
            {'from': 'aspirin', 'to': 'ibuprofen', 'relationship': 'SIMILAR_TO', 'weight': 0.7},
            {'from': 'res_1', 'to': 'aspirin', 'relationship': 'CONTRIBUTED_TO', 'weight': 1.0},
            {'from': 'res_2', 'to': 'caffeine', 'relationship': 'CONTRIBUTED_TO', 'weight': 1.0},
            {'from': 'res_3', 'to': 'vitamin_c', 'relationship': 'CONTRIBUTED_TO', 'weight': 1.0},
        ])
    
    def query(self, query_type, **kwargs):
        """Query the knowledge graph"""
        if query_type == 'find_similar':
            molecule_id = kwargs.get('molecule_id')
            similar = []
            for edge in self.edges:
                if edge['relationship'] == 'SIMILAR_TO' and edge['from'] == molecule_id:
                    similar.append({
                        'molecule_id': edge['to'],
                        'similarity': edge['weight']
                    })
            return similar
        
        elif query_type == 'researcher_network':
            researcher_id = kwargs.get('researcher_id')
            contributions = []
            for edge in self.edges:
                if edge['relationship'] == 'CONTRIBUTED_TO' and edge['from'] == researcher_id:
                    contributions.append({
                        'molecule_id': edge['to'],
                        'contribution_weight': edge['weight']
                    })
            return contributions
        
        elif query_type == 'molecule_properties':
            molecule_id = kwargs.get('molecule_id')
            properties = []
            for edge in self.edges:
                if edge['relationship'] == 'HAS_PROPERTY' and edge['from'] == molecule_id:
                    properties.append({
                        'property_id': edge['to'],
                        'confidence': edge['weight']
                    })
            return properties
        
        return []

knowledge_graph = KnowledgeGraph()

# ===== WEB SCRAPING AGENT =====
class WebScrapingAgent:
    def __init__(self):
        self.active_jobs = []
        self.results_cache = {}
        
    def search_pubchem(self, query):
        """Simulate PubChem API search"""
        time.sleep(0.5)  # Simulate API delay
        
        # Mock response
        return {
            'source': 'PubChem',
            'query': query,
            'results': [
                {
                    'cid': random.randint(1000, 9999),
                    'name': query.title(),
                    'formula': 'C9H8O4' if 'aspirin' in query.lower() else 'C8H10N4O2',
                    'molecular_weight': random.randint(100, 500),
                    'synonyms': [f"{query} synonym {i}" for i in range(3)],
                    'url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{random.randint(1000, 9999)}"
                }
            ],
            'timestamp': datetime.now().isoformat() + "Z"
        }
    
    def search_pubmed(self, query):
        """Simulate PubMed API search"""
        time.sleep(0.7)
        
        return {
            'source': 'PubMed',
            'query': query,
            'results': [
                {
                    'pmid': f"PMC{random.randint(1000000, 9999999)}",
                    'title': f"Study of {query} in clinical trials",
                    'authors': ['Smith, J.', 'Chen, R.', 'Garcia, M.'],
                    'journal': 'Nature Medicine',
                    'year': random.randint(2018, 2024),
                    'abstract': f"This study investigates the effects of {query} on various biological pathways...",
                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{random.randint(30000000, 39999999)}/"
                }
                for _ in range(random.randint(1, 3))
            ],
            'timestamp': datetime.now().isoformat() + "Z"
        }
    
    def search_clinical_trials(self, query):
        """Simulate clinical trials search"""
        time.sleep(0.9)
        
        return {
            'source': 'ClinicalTrials.gov',
            'query': query,
            'results': [
                {
                    'nct_id': f"NCT{random.randint(10000000, 99999999)}",
                    'title': f"Phase {random.randint(1, 3)} Trial of {query}",
                    'status': random.choice(['Recruiting', 'Completed', 'Active']),
                    'conditions': ['Condition A', 'Condition B'],
                    'intervention': query,
                    'url': f"https://clinicaltrials.gov/ct2/show/NCT{random.randint(10000000, 99999999)}"
                }
            ],
            'timestamp': datetime.now().isoformat() + "Z"
        }
    
    def search_all(self, query):
        """Search all sources concurrently"""
        print(f"🔍 Web agent searching for: {query}")
        
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {
                executor.submit(self.search_pubchem, query): 'pubchem',
                executor.submit(self.search_pubmed, query): 'pubmed',
                executor.submit(self.search_clinical_trials, query): 'clinical_trials'
            }
            
            results = {}
            for future in as_completed(futures):
                source = futures[future]
                try:
                    results[source] = future.result(timeout=10)
                except Exception as e:
                    results[source] = {'error': str(e), 'source': source}
            
            # Create scraping job record
            job_id = f"job_{len(web_scraping_jobs) + 1:04d}"
            job = {
                'id': job_id,
                'query': query,
                'status': 'completed',
                'results': results,
                'timestamp': datetime.now().isoformat() + "Z"
            }
            
            web_scraping_jobs.append(job)
            self.results_cache[query] = job
            
            return job

web_agent = WebScrapingAgent()

# ===== INITIAL SAMPLE DATA =====
def init_sample_data():
    global molecules_db, contributions_db, researchers_db
    
    researchers_db.extend([
        {
            "id": "res_001",
            "name": "Dr. Jane Smith",
            "email": "jane.smith@oxford.edu",
            "institution": "University of Oxford",
            "orcid": "0000-0001-2345-6789",
            "expertise": ["Computational Chemistry", "Drug Discovery"],
            "reputation_score": 95,
            "molecules_contributed": 12,
            "total_contributions": 47,
            "wallet_address": "0x742d35Cc6634C0532925a3b844Bc9e",
            "joined": "2023-06-15T09:30:00Z"
        },
        {
            "id": "res_002",
            "name": "Dr. Robert Chen",
            "email": "robert.chen@cambridge.edu",
            "institution": "University of Cambridge",
            "orcid": "0000-0002-3456-7890",
            "expertise": ["Medicinal Chemistry", "AI/ML"],
            "reputation_score": 88,
            "molecules_contributed": 8,
            "total_contributions": 32,
            "wallet_address": "0x2aC3C5C5F5F5F5F5F5F5F5F5F5F5F5F5",
            "joined": "2023-08-22T14:20:00Z"
        },
        {
            "id": "res_003",
            "name": "Dr. Maria Garcia",
            "email": "maria.garcia@ethz.ch",
            "institution": "ETH Zurich",
            "orcid": "0000-0003-4567-8901",
            "expertise": ["Structural Biology", "Biochemistry"],
            "reputation_score": 92,
            "molecules_contributed": 15,
            "total_contributions": 56,
            "wallet_address": "0x3bD4E6F8A0B1C2D3E4F5A6B7C8D9E0F",
            "joined": "2023-05-10T11:15:00Z"
        }
    ])
    
    molecules_db.extend([
        {
            "id": "mol_001",
            "name": "Aspirin",
            "formula": "C9H8O4",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "molecular_weight": 180.16,
            "description": "Nonsteroidal anti-inflammatory drug used for pain relief, fever reduction, and anti-inflammatory purposes.",
            "ai_score": 87,
            "consensus_breakdown": {
                "toxicity": 0.92,
                "solubility": 0.78,
                "binding": 0.81,
                "druglikeness": 0.88,
                "synthesis": 0.76
            },
            "on_blockchain": True,
            "blockchain_data": {
                "tx_hash": "0x7a3f9c8d2e1b4a5f6c7d8e9f0a1b2c3d4e5f67890",
                "block_number": 1245678,
                "value_usd": 1245.67,
                "registered_by": "res_001",
                "registration_date": "2024-01-15T10:30:00Z"
            },
            "properties": {
                "toxicity": {"score": 0.92, "category": "Low"},
                "solubility": {"score": 0.78, "category": "Moderate", "logS": -1.5},
                "melting_point": {"value": 135, "unit": "°C"},
                "bioavailability": {"score": 0.85, "category": "High"},
                "logP": 1.19,
                "polar_surface_area": 63.6,
                "rotatable_bonds": 3,
                "h_bond_donors": 1,
                "h_bond_acceptors": 4
            },
            "applications": ["Pain relief", "Anti-inflammatory", "Antiplatelet"],
            "contribution_count": 3,
            "created_at": "2024-01-15T10:30:00Z",
            "created_by": "res_001",
            "tags": ["NSAID", "Analgesic", "Common Drug"],
            "3d_structure_url": "/api/molecules/mol_001/structure",
            "similar_molecules": ["ibuprofen", "naproxen"]
        },
        {
            "id": "mol_002",
            "name": "Caffeine",
            "formula": "C8H10N4O2",
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
            "molecular_weight": 194.19,
            "description": "Central nervous system stimulant of the methylxanthine class.",
            "ai_score": 92,
            "consensus_breakdown": {
                "toxicity": 0.85,
                "solubility": 0.82,
                "binding": 0.88,
                "druglikeness": 0.94,
                "synthesis": 0.89
            },
            "on_blockchain": True,
            "blockchain_data": {
                "tx_hash": "0x8b4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e",
                "block_number": 1245679,
                "value_usd": 892.45,
                "registered_by": "res_002",
                "registration_date": "2024-01-16T14:20:00Z"
            },
            "properties": {
                "toxicity": {"score": 0.85, "category": "Moderate"},
                "solubility": {"score": 0.82, "category": "Moderate", "logS": -0.8},
                "melting_point": {"value": 238, "unit": "°C"},
                "bioavailability": {"score": 0.91, "category": "High"},
                "logP": -0.07,
                "polar_surface_area": 58.4,
                "rotatable_bonds": 0,
                "h_bond_donors": 0,
                "h_bond_acceptors": 6
            },
            "applications": ["Stimulant", "Alertness", "Performance Enhancer"],
            "contribution_count": 5,
            "created_at": "2024-01-16T14:20:00Z",
            "created_by": "res_002",
            "tags": ["Stimulant", "Alkaloid", "Common"],
            "3d_structure_url": "/api/molecules/mol_002/structure",
            "similar_molecules": ["theophylline", "theobromine"]
        }
    ])
    
    contributions_db.extend([
        {
            "id": "cont_001",
            "molecule_id": "mol_001",
            "type": "quantum_analysis",
            "title": "DFT Study of Aspirin's Electronic Properties",
            "researcher_id": "res_001",
            "description": "Density functional theory calculations reveal novel electronic properties of aspirin.",
            "methodology": "B3LYP/6-311+G(d,p) level of theory with implicit solvent model",
            "results": "Identified new charge transfer pathways and binding affinities",
            "data_url": "https://example.com/data/aspirin_dft.zip",
            "ai_validation": 94,
            "files": ["optimized_geometry.xyz", "electronic_spectra.csv", "binding_energies.json"],
            "timestamp": "2024-01-20T11:30:00Z",
            "blockchain_tx": "0x9c5f6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c",
            "reputation_reward": 25,
            "status": "validated"
        },
        {
            "id": "cont_002",
            "molecule_id": "mol_002",
            "type": "experimental",
            "title": "Solubility Enhancement of Caffeine in Various Solvents",
            "researcher_id": "res_002",
            "description": "Experimental investigation of caffeine solubility enhancement techniques.",
            "methodology": "UV-Vis spectroscopy, HPLC analysis, particle size reduction",
            "results": "45% solubility improvement using co-solvent approach",
            "data_url": "https://example.com/data/caffeine_solubility.zip",
            "ai_validation": 88,
            "files": ["spectra_data.csv", "hpc_chromatograms.pdf", "particle_size_analysis.json"],
            "timestamp": "2024-01-19T14:20:00Z",
            "blockchain_tx": "0xa1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0",
            "reputation_reward": 18,
            "status": "validated"
        }
    ])
    
    # Initialize blockchain transactions
    for molecule in molecules_db:
        if molecule.get('on_blockchain') and molecule.get('blockchain_data'):
            blockchain_txs.append(molecule['blockchain_data'])

init_sample_data()

# ===== REPORT GENERATION FUNCTIONS =====

def generate_analysis_report(molecule_data, analysis_results):
    """Generate a comprehensive AI analysis report"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    report = f"""
╔════════════════════════════════════════════════════════════════════════════╗
║                    AI ANALYSIS REPORT - MOLCHAIN                          ║
╚════════════════════════════════════════════════════════════════════════════╝

Generated: {timestamp}
Report ID: {molecule_data.get('id', 'UNKNOWN')}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

📋 MOLECULE INFORMATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Name:                {molecule_data.get('name', 'N/A')}
  SMILES:              {molecule_data.get('smiles', 'N/A')}
  Molecular Formula:   {molecule_data.get('formula', 'N/A')}
  Molecular Weight:    {molecule_data.get('molecular_weight', 'N/A')} g/mol
  Description:         {molecule_data.get('description', 'No description provided')}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

🤖 AI AGENT ANALYSIS RESULTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

"""
    
    if 'agent_results' in analysis_results:
        for agent_id, result in analysis_results['agent_results'].items():
            agent_name = result.get('agent_name', agent_id)
            score = result.get('score', 0.0)
            confidence = result.get('confidence', 0.0)
            
            # Create visual bar
            bar_length = int(score * 20)
            bar = '█' * bar_length + '░' * (20 - bar_length)
            
            report += f"""
  {agent_name}
    Score:      {score} [{bar}] {int(score*100)}%
    Confidence: {confidence}
    Status:     {'✓ Success' if result.get('status') == 'success' else '✗ Failed'}
"""
    
    report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

📊 CONSENSUS ANALYSIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Consensus Score:     {analysis_results.get('consensus_score', 0.0)}
  Consensus %:         {int(analysis_results.get('consensus_score', 0.0)*100)}%
  
  Analysis Status:     {'✓ PASSED' if analysis_results.get('blockchain_ready') else '✗ FAILED'}
  Blockchain Ready:    {'Yes' if analysis_results.get('blockchain_ready') else 'No'}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

🔗 BLOCKCHAIN REGISTRATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Network:             Flare Coston (Chain ID: 16)
  Status:              {'Registered' if molecule_data.get('on_blockchain') else 'Pending'}
"""
    
    if molecule_data.get('blockchain_data'):
        blockchain_data = molecule_data.get('blockchain_data', {})
        report += f"""
  TX Hash:             {blockchain_data.get('tx_hash', 'N/A')}
  Block Number:        {blockchain_data.get('block_number', 'N/A')}
  Gas Used:            {blockchain_data.get('gas_used', 'N/A')} gwei
  Confirmation:        ✓ Confirmed
  Value (USD):         ${blockchain_data.get('value_usd', 0.0)}
  Explorer URL:        {blockchain_data.get('explorer_url', 'N/A')}
"""
    
    report += """
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

📈 RECOMMENDATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

"""
    
    consensus_score = analysis_results.get('consensus_score', 0.0)
    
    if consensus_score > 0.8:
        report += """
  1. ✓ Excellent drug candidate - proceed with synthesis planning
  2. ✓ High confidence in all AI agent predictions
  3. ✓ Recommended for blockchain registration and IP protection
  4. ✓ Consider for clinical trial pipeline
  5. ✓ Strong market potential identified
"""
    elif consensus_score > 0.6:
        report += """
  1. ✓ Good drug candidate - further optimization recommended
  2. ✓ Moderate confidence in predictions
  3. ✓ Recommended for additional computational analysis
  4. ✓ Consider for further binding affinity studies
  5. ⚠ May require solubility improvements
"""
    else:
        report += """
  1. ⚠ Moderate drug candidate - significant refinement needed
  2. ⚠ Lower confidence in current predictions
  3. ✓ Recommended for structural modifications
  4. ⚠ Consider alternative scaffolds
  5. ⚠ Further analysis required before advancement
"""
    
    report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

⚙️ ANALYSIS METADATA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Analysis Timestamp:  {analysis_results.get('timestamp', 'N/A')}
  AI Agents Count:     5
  Processing Status:   ✓ Complete
  Report Version:      1.0

╔════════════════════════════════════════════════════════════════════════════╗
║  This report was generated by the MolChain AI System powered by Flare      ║
║  Blockchain. All analysis results are simulated by AI agents with 24/7     ║
║  consensus algorithms for maximum accuracy and reliability.                ║
╚════════════════════════════════════════════════════════════════════════════╝
"""
    
    return report


# ===== ENHANCED API ENDPOINTS =====

@app.route('/api/health', methods=['GET'])
@rate_limit(max_per_minute=120)
def health():
    """Enhanced health check with swarm status"""
    return jsonify({
        "status": "healthy",
        "service": "MolChain AI + Blockchain Backend",
        "version": app.config['VERSION'],
        "timestamp": datetime.now().isoformat() + "Z",
        "stats": {
            "molecules": len(molecules_db),
            "contributions": len(contributions_db),
            "researchers": len(researchers_db),
            "ai_agents": len(ai_agents_db),
            "blockchain_txs": len(blockchain_txs),
            "knowledge_graph_nodes": len(knowledge_graph.nodes),
            "web_scraping_jobs": len(web_scraping_jobs)
        },
        "capabilities": [
            "AI Agent Swarm (5 agents)",
            "Flare Blockchain Simulation",
            "Knowledge Graph",
            "Web Scraping Agent",
            "Real-time Analytics",
            "Blockchain Integration"
        ],
        "uptime": "100%"
    })

@app.route('/api/molecules', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_molecules():
    """Enhanced molecule listing with filters"""
    try:
        page = int(request.args.get('page', 1))
        limit = int(request.args.get('limit', 20))
        sort_by = request.args.get('sort', 'ai_score')
        filter_type = request.args.get('filter', 'all')
        search = request.args.get('search', '').lower()
        
        # Filter molecules
        filtered = molecules_db.copy()
        
        if filter_type == 'verified':
            filtered = [m for m in filtered if m.get('on_blockchain')]
        elif filter_type == 'pending':
            filtered = [m for m in filtered if not m.get('on_blockchain')]
        
        # Search
        if search:
            filtered = [m for m in filtered if 
                       search in m.get('name', '').lower() or
                       search in m.get('formula', '').lower() or
                       search in m.get('smiles', '').lower() or
                       any(search in tag.lower() for tag in m.get('tags', []))]
        
        # Sort
        reverse = True
        if sort_by == 'ai_score':
            filtered.sort(key=lambda x: x.get('ai_score', 0), reverse=True)
        elif sort_by == 'name':
            filtered.sort(key=lambda x: x.get('name', ''))
            reverse = False
        elif sort_by == 'weight':
            filtered.sort(key=lambda x: x.get('molecular_weight', 0))
        elif sort_by == 'recent':
            filtered.sort(key=lambda x: x.get('created_at', ''), reverse=True)
        elif sort_by == 'value':
            filtered.sort(key=lambda x: x.get('blockchain_data', {}).get('value_usd', 0), reverse=True)
        
        # Paginate
        start_idx = (page - 1) * limit
        end_idx = start_idx + limit
        paginated = filtered[start_idx:end_idx]
        
        return jsonify({
            "success": True,
            "count": len(filtered),
            "page": page,
            "limit": limit,
            "total_pages": (len(filtered) + limit - 1) // limit,
            "molecules": paginated,
            "filters_applied": {
                "sort": sort_by,
                "filter": filter_type,
                "search": search if search else None
            }
        })
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/molecules/<molecule_id>', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_molecule(molecule_id):
    """Enhanced single molecule with all related data"""
    molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
    if not molecule:
        return jsonify({"success": False, "error": "Molecule not found"}), 404
    
    # Get related data
    molecule_contributions = [c for c in contributions_db if c["molecule_id"] == molecule_id]
    similar_molecules = knowledge_graph.query('find_similar', molecule_id=molecule_id)
    
    # Get agent status
    agent_status = []
    for agent in ai_agents_db:
        if agent['status'] == 'active':
            agent_status.append({
                "name": agent['name'],
                "status": "active",
                "accuracy": agent['accuracy'],
                "tasks_completed": agent['tasks_completed'],
                "last_active": agent['last_active']
            })
    
    # Get knowledge graph connections
    kg_connections = knowledge_graph.query('molecule_properties', molecule_id=molecule_id)
    
    # Get molecule value from blockchain
    molecule_value = flare_simulator.get_molecule_value(molecule_id, molecule.get('ai_score', 0))
    
    return jsonify({
        "success": True,
        "molecule": molecule,
        "contributions": molecule_contributions,
        "ai_agents": agent_status,
        "knowledge_graph": {
            "connections": kg_connections,
            "similar_molecules": similar_molecules
        },
        "blockchain": {
            "value_usd": molecule_value,
            "ftso_price": flare_simulator.ftso_price,
            "on_chain": molecule.get('on_blockchain', False),
            "tx_hash": molecule.get('blockchain_data', {}).get('tx_hash')
        },
        "analytics": {
            "total_views": random.randint(100, 1000),
            "citation_count": random.randint(5, 50),
            "research_interest": "High" if molecule.get('ai_score', 0) > 80 else "Moderate"
        }
    })

@app.route('/api/molecules', methods=['POST'])
@rate_limit(max_per_minute=30)
def create_molecule():
    """Enhanced molecule creation with AI validation"""
    try:
        data = request.json
        if not data.get('name') or not data.get('smiles'):
            return jsonify({"success": False, "error": "Missing required fields: name and smiles"}), 400
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(data['smiles'])
        if not mol:
            return jsonify({"success": False, "error": "Invalid SMILES string"}), 400
        
        # Generate molecule ID
        molecule_id = f"mol_{len(molecules_db) + 1:03d}"
        
        # Calculate basic properties
        mol_wt = Descriptors.MolLogP(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        
        new_molecule = {
            "id": molecule_id,
            "name": data['name'],
            "formula": formula,
            "smiles": data['smiles'],
            "inchi": Chem.MolToInchi(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "description": data.get('description', ''),
            "ai_score": 0,
            "consensus_breakdown": None,
            "on_blockchain": False,
            "blockchain_data": None,
            "properties": {
                "toxicity": {"score": 0.5, "category": "Pending"},
                "solubility": {"score": 0.5, "category": "Pending"},
                "melting_point": {"value": None, "unit": "°C"},
                "bioavailability": {"score": 0.5, "category": "Pending"},
                "logP": round(Descriptors.MolLogP(mol), 2),
                "polar_surface_area": round(Descriptors.TPSA(mol), 1),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "h_bond_donors": Descriptors.NumHDonors(mol),
                "h_bond_acceptors": Descriptors.NumHAcceptors(mol)
            },
            "applications": data.get('applications', []),
            "contribution_count": 0,
            "created_at": datetime.now().isoformat() + "Z",
            "created_by": data.get('researcher_id', 'anonymous'),
            "tags": data.get('tags', []),
            "3d_structure_url": f"/api/molecules/{molecule_id}/structure",
            "similar_molecules": []
        }
        
        molecules_db.append(new_molecule)
        
        # Start AI analysis in background
        def analyze_in_background():
            analysis = ai_swarm.analyze_molecule(data['smiles'], molecule_id)
            new_molecule['ai_score'] = round(analysis['consensus_score'] * 100)
            new_molecule['consensus_breakdown'] = {
                agent_id: result['score']
                for agent_id, result in analysis['agent_results'].items()
            }
            
            # Update properties based on AI results
            if 'toxicity' in analysis['agent_results']:
                tox_score = analysis['agent_results']['toxicity']['score']
                new_molecule['properties']['toxicity'] = {
                    "score": tox_score,
                    "category": "Low" if tox_score > 0.8 else "Moderate" if tox_score > 0.5 else "High"
                }
            
            if 'solubility' in analysis['agent_results']:
                sol_score = analysis['agent_results']['solubility']['score']
                new_molecule['properties']['solubility'] = {
                    "score": sol_score,
                    "category": "High" if sol_score > 0.8 else "Moderate" if sol_score > 0.5 else "Low",
                    "logS": round(-2 * (1 - sol_score), 2)
                }
        
        # Run analysis in background thread
        threading.Thread(target=analyze_in_background).start()
        
        return jsonify({
            "success": True,
            "message": "Molecule created successfully. AI analysis started.",
            "molecule": new_molecule,
            "analysis_in_progress": True,
            "estimated_completion": "30 seconds"
        }), 201
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

# ===== HELPER FUNCTIONS =====
def generate_recommendations(analysis):
    """Generate recommendations based on AI analysis"""
    recommendations = []
    
    # Check toxicity
    if 'toxicity' in analysis['agent_results']:
        tox_score = analysis['agent_results']['toxicity']['score']
        if tox_score < 0.3:
            recommendations.append("⚠️ High toxicity detected - consider structural modifications")
        elif tox_score < 0.6:
            recommendations.append("🔬 Moderate toxicity - further testing recommended")
    
    # Check solubility
    if 'solubility' in analysis['agent_results']:
        sol_score = analysis['agent_results']['solubility']['score']
        if sol_score < 0.4:
            recommendations.append("💧 Low solubility - consider formulation improvements")
    
    # Check drug likeness
    if 'druglikeness' in analysis['agent_results']:
        dl_score = analysis['agent_results']['druglikeness']['score']
        if dl_score < 0.5:
            recommendations.append("💊 Poor drug-likeness - review molecular properties")
    
    return recommendations

@app.route('/api/analyze', methods=['POST'])
@rate_limit(max_per_minute=20)
def analyze_molecule():
    """Enhanced AI analysis endpoint"""
    try:
        data = request.json
        molecule_id = data.get('molecule_id')
        smiles = data.get('smiles')
        
        if not molecule_id and not smiles:
            return jsonify({"success": False, "error": "Provide either molecule_id or smiles"}), 400
        
        # Get molecule if ID provided
        molecule = None
        if molecule_id:
            molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
            if not molecule:
                return jsonify({"success": False, "error": "Molecule not found"}), 404
            smiles = molecule['smiles']
        
        # Run AI analysis
        print(f"🔬 Starting comprehensive AI analysis for {smiles}")
        analysis = ai_swarm.analyze_molecule(smiles, molecule_id or 'new_molecule')
        
        # Register on blockchain if consensus is high
        blockchain_data = None
        if analysis['blockchain_ready'] and molecule:
            blockchain_data = flare_simulator.register_molecule(molecule, analysis)
            molecule['on_blockchain'] = True
            molecule['blockchain_data'] = blockchain_data
        
        response = {
            "success": True,
            "message": "AI analysis completed successfully",
            "consensus_score": analysis['consensus_score'],
            "blockchain_ready": analysis['blockchain_ready'],
            "agent_results": analysis['agent_results'],
            "blockchain": blockchain_data,
            "recommendations": generate_recommendations(analysis),
            "timestamp": analysis['timestamp']
        }
        
        if molecule:
            response['molecule'] = {
                'id': molecule['id'],
                'name': molecule['name'],
                'ai_score': round(analysis['consensus_score'] * 100),
                'on_blockchain': molecule['on_blockchain']
            }
        
        return jsonify(response)
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/analyze/report', methods=['POST'])
@rate_limit(max_per_minute=20)
def download_analysis_report():
    """Generate and download AI analysis report"""
    try:
        data = request.json
        molecule_id = data.get('molecule_id')
        
        if not molecule_id:
            return jsonify({"success": False, "error": "molecule_id is required"}), 400
        
        # Get molecule from database
        molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
        if not molecule:
            return jsonify({"success": False, "error": "Molecule not found"}), 404
        
        # Get analysis (run if not already done)
        smiles = molecule.get('smiles')
        analysis = ai_swarm.analyze_molecule(smiles, molecule_id)
        
        # Update molecule with blockchain data if ready
        if analysis['blockchain_ready'] and not molecule.get('on_blockchain'):
            blockchain_data = flare_simulator.register_molecule(molecule, analysis)
            molecule['on_blockchain'] = True
            molecule['blockchain_data'] = blockchain_data
        
        # Generate report
        report = generate_analysis_report(molecule, analysis)
        
        # Create downloadable file
        output = BytesIO()
        output.write(report.encode('utf-8'))
        output.seek(0)
        
        # Generate filename
        filename = f"AI_Analysis_Report_{molecule_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        
        return send_file(
            output,
            mimetype='text/plain',
            as_attachment=True,
            download_name=filename
        )
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/analyze/batch', methods=['POST'])
@rate_limit(max_per_minute=10)
def analyze_batch():
    """Analyze multiple molecules in batch"""
    try:
        data = request.json
        smiles_list = data.get('smiles', [])
        
        if not smiles_list or len(smiles_list) > 10:
            return jsonify({
                "success": False, 
                "error": "Provide up to 10 SMILES strings"
            }), 400
        
        results = []
        for smiles in smiles_list:
            analysis = ai_swarm.analyze_molecule(smiles, f"batch_{hash(smiles)}")
            results.append({
                "smiles": smiles,
                "consensus_score": analysis['consensus_score'],
                "blockchain_ready": analysis['blockchain_ready'],
                "top_agents": {
                    agent_id: result['score']
                    for agent_id, result in list(analysis['agent_results'].items())[:3]
                }
            })
        
        return jsonify({
            "success": True,
            "count": len(results),
            "results": results,
            "batch_id": f"batch_{int(time.time())}",
            "timestamp": datetime.now().isoformat() + "Z"
        })
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/stats', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_stats():
    """Enhanced statistics endpoint"""
    total_ai_score = sum(m.get('ai_score', 0) for m in molecules_db)
    avg_ai_score = round(total_ai_score / len(molecules_db)) if molecules_db else 0
    
    # Calculate agent statistics
    agent_stats = []
    for agent in ai_agents_db:
        agent_stats.append({
            "name": agent['name'],
            "tasks": agent['tasks_completed'],
            "accuracy": agent['accuracy'],
            "status": agent['status']
        })
    
    # Calculate blockchain value
    total_value = sum(
        m.get('blockchain_data', {}).get('value_usd', 0) 
        for m in molecules_db 
        if m.get('on_blockchain')
    )
    
    # Get recent activity
    recent_molecules = sorted(
        molecules_db, 
        key=lambda x: x.get('created_at', ''), 
        reverse=True
    )[:5]
    
    recent_contributions = sorted(
        contributions_db,
        key=lambda x: x.get('timestamp', ''),
        reverse=True
    )[:5]
    
    return jsonify({
        "success": True,
        "stats": {
            "total_molecules": len(molecules_db),
            "verified_molecules": len([m for m in molecules_db if m.get('on_blockchain')]),
            "total_contributions": len(contributions_db),
            "total_researchers": len(researchers_db),
            "avg_ai_score": avg_ai_score,
            "ai_agents_active": len([a for a in ai_agents_db if a['status'] == 'active']),
            "total_ai_analyses": sum(a['tasks_completed'] for a in ai_agents_db),
            "blockchain_txs": len(blockchain_txs),
            "total_value_locked": round(total_value, 2),
            "knowledge_graph_nodes": len(knowledge_graph.nodes),
            "knowledge_graph_edges": len(knowledge_graph.edges),
            "web_scraping_jobs": len(web_scraping_jobs)
        },
        "agent_statistics": agent_stats,
        "recent_activity": {
            "molecules": [{"id": m['id'], "name": m['name'], "created": m.get('created_at')} 
                         for m in recent_molecules],
            "contributions": [{"id": c['id'], "title": c['title'], "timestamp": c.get('timestamp')} 
                            for c in recent_contributions]
        },
        "performance": {
            "avg_response_time": "0.45s",
            "uptime": "99.9%",
            "concurrent_analyses": app.config['AI_AGENTS'],
            "cache_hit_rate": "92%"
        },
        "timestamp": datetime.now().isoformat() + "Z"
    })

@app.route('/api/agents', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_agents():
    """Enhanced AI agents endpoint"""
    return jsonify({
        "success": True,
        "swarm": {
            "name": "MolChain AI Agent Swarm v2.0",
            "status": "active",
            "version": "2.0.0-hackathon",
            "architecture": "Multi-agent consensus system",
            "total_agents": app.config['AI_AGENTS'],
            "active_agents": len([a for a in ai_agents_db if a['status'] == 'active']),
            "total_tasks": sum(a['tasks_completed'] for a in ai_agents_db),
            "consensus_accuracy": 0.89,
            "avg_response_time": "1.2s",
            "last_consensus_check": ai_swarm.consensus_history[-1]['timestamp'] if ai_swarm.consensus_history else None
        },
        "agents": ai_agents_db,
        "consensus_history": ai_swarm.consensus_history[-10:],  # Last 10 consensus results
        "capabilities": [
            "Parallel multi-agent analysis",
            "Weighted consensus scoring",
            "Real-time blockchain validation",
            "Cross-validation with knowledge graph",
            "Automated recommendation generation"
        ]
    })

@app.route('/api/agents/<agent_id>/status', methods=['GET'])
def get_agent_status(agent_id):
    """Get detailed status of a specific AI agent"""
    agent = next((a for a in ai_agents_db if a['id'] == agent_id), None)
    if not agent:
        return jsonify({"success": False, "error": "Agent not found"}), 404
    
    # Get recent tasks for this agent
    recent_tasks = []
    for consensus in ai_swarm.consensus_history[-20:]:
        if agent_id in consensus.get('agent_results', {}):
            task_result = consensus['agent_results'][agent_id]
            recent_tasks.append({
                "molecule_id": consensus['molecule_id'],
                "smiles": consensus['smiles'],
                "score": task_result.get('score'),
                "confidence": task_result.get('confidence'),
                "timestamp": consensus['timestamp']
            })
    
    return jsonify({
        "success": True,
        "agent": agent,
        "performance": {
            "tasks_last_hour": len([t for t in recent_tasks if 
                                   datetime.fromisoformat(t['timestamp'].replace('Z', '')) > 
                                   datetime.now() - timedelta(hours=1)]),
            "avg_score": round(sum(t.get('score', 0) for t in recent_tasks) / len(recent_tasks), 3) if recent_tasks else 0,
            "success_rate": 0.98,
            "avg_processing_time": "0.8s"
        },
        "recent_tasks": recent_tasks[:10],
        "resource_usage": {
            "memory_mb": random.randint(50, 200),
            "cpu_percent": random.randint(1, 30),
            "gpu_available": torch.cuda.is_available(),
            "model_size_mb": random.randint(100, 500)
        }
    })

@app.route('/api/flare/price', methods=['GET'])
@rate_limit(max_per_minute=120)
def get_flare_price():
    """Enhanced Flare FTSO price feed"""
    price_data = flare_simulator.get_ftso_price()
    
    return jsonify({
        "success": True,
        "oracle": "Flare Time Series Oracle (FTSO)",
        "network": "Coston Testnet",
        "chain_id": flare_simulator.chain_id,
        "data": price_data,
        "endpoints": {
            "ftso_manager": "0x1000000000000000000000000000000000000003",
            "ftso_reward_manager": "0x1000000000000000000000000000000000000004",
            "price_submitter": "0x1000000000000000000000000000000000000005"
        },
        "protocols_used": ["FTSO", "State Connector", "FAssets"],
        "update_frequency": "Every 90 seconds",
        "decimals": 18,
        "confidence_threshold": 0.85
    })

@app.route('/api/flare/history', methods=['GET'])
def get_flare_price_history():
    """Get Flare price history"""
    hours = int(request.args.get('hours', 24))
    data_points = min(hours * 4, 100)  # Max 100 points
    
    history = []
    base_time = datetime.now()
    base_price = flare_simulator.ftso_price
    
    for i in range(data_points):
        time_offset = timedelta(minutes=-15 * i)
        price_variation = random.uniform(-0.002, 0.002)
        price = max(0.01, base_price + price_variation)
        
        history.append({
            "timestamp": (base_time + time_offset).isoformat() + "Z",
            "price_usd": round(price, 4),
            "confidence": round(random.uniform(0.8, 0.99), 3),
            "volume": random.randint(1000000, 5000000)
        })
    
    return jsonify({
        "success": True,
        "timeframe_hours": hours,
        "data_points": len(history),
        "current_price": round(flare_simulator.ftso_price, 4),
        "24h_change": round(random.uniform(-0.05, 0.05), 3),
        "history": history
    })

@app.route('/api/flare/<molecule_id>', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_molecule_value(molecule_id):
    """Enhanced molecular valuation from Flare oracle"""
    molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
    if not molecule:
        return jsonify({"success": False, "error": "Molecule not found"}), 404
    
    ai_score = molecule.get('ai_score', 0) / 100.0
    value = flare_simulator.get_molecule_value(molecule_id, ai_score)
    
    # Get price data
    price_data = flare_simulator.get_ftso_price()
    
    return jsonify({
        "success": True,
        "molecule": molecule_id,
        "name": molecule.get('name'),
        "valuation": {
            "value_usd": value,
            "valuation_model": "AI Consensus × Market Factor × FTSO Price",
            "components": {
                "ai_consensus_score": ai_score,
                "market_factor": round(random.uniform(0.8, 1.2), 3),
                "ftso_price_usd": price_data['price_usd'],
                "molecular_complexity": random.uniform(0.5, 2.0)
            },
            "breakdown": {
                "base_value": round(ai_score * 1000, 2),
                "market_adjustment": round(random.uniform(0.8, 1.2), 3),
                "oracle_price": price_data['price_usd'],
                "final_value": value
            }
        },
        "oracle": {
            "name": "Flare FTSO Enhanced",
            "confidence": price_data['confidence'],
            "timestamp": price_data['timestamp'],
            "update_interval": "90 seconds"
        },
        "blockchain": {
            "on_chain": molecule.get('on_blockchain', False),
            "tx_hash": molecule.get('blockchain_data', {}).get('tx_hash'),
            "block_number": molecule.get('blockchain_data', {}).get('block_number'),
            "registered_value": molecule.get('blockchain_data', {}).get('value_usd')
        } if molecule.get('on_blockchain') else None,
        "timestamp": datetime.now().isoformat() + "Z"
    })

@app.route('/api/flare/register', methods=['POST'])
@rate_limit(max_per_minute=20)
def register_on_blockchain():
    """Register molecule on Flare blockchain"""
    try:
        data = request.json
        molecule_id = data.get('molecule_id')
        
        if not molecule_id:
            return jsonify({"success": False, "error": "Missing molecule_id"}), 400
        
        molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
        if not molecule:
            return jsonify({"success": False, "error": "Molecule not found"}), 404
        
        if molecule.get('on_blockchain'):
            return jsonify({
                "success": False, 
                "error": "Molecule already registered on blockchain"
            }), 400
        
        # Ensure molecule has AI score
        if molecule.get('ai_score', 0) == 0:
            return jsonify({
                "success": False,
                "error": "Molecule needs AI analysis before blockchain registration"
            }), 400
        
        # Register on blockchain
        ai_scores = {
            'consensus_score': molecule.get('ai_score', 0) / 100.0
        }
        
        tx_data = flare_simulator.register_molecule(molecule, ai_scores)
        
        # Update molecule
        molecule['on_blockchain'] = True
        molecule['blockchain_data'] = tx_data
        
        return jsonify({
            "success": True,
            "message": "Molecule successfully registered on Flare blockchain",
            "transaction": tx_data,
            "molecule": {
                "id": molecule_id,
                "name": molecule.get('name'),
                "value_usd": tx_data['value_usd']
            },
            "next_steps": [
                "View transaction on explorer",
                "Share with research community",
                "Monitor valuation changes"
            ]
        })
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/contributions', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_contributions():
    """Enhanced contributions endpoint"""
    molecule_id = request.args.get('molecule_id')
    researcher_id = request.args.get('researcher_id')
    contribution_type = request.args.get('type')
    
    filtered = contributions_db.copy()
    
    if molecule_id:
        filtered = [c for c in filtered if c["molecule_id"] == molecule_id]
    if researcher_id:
        filtered = [c for c in filtered if c["researcher_id"] == researcher_id]
    if contribution_type:
        filtered = [c for c in filtered if c["type"] == contribution_type]
    
    # Sort by AI validation score
    filtered.sort(key=lambda x: x.get('ai_validation', 0), reverse=True)
    
    return jsonify({
        "success": True,
        "count": len(filtered),
        "contributions": filtered,
        "summary": {
            "types": list(set(c['type'] for c in filtered)),
            "avg_validation": round(sum(c.get('ai_validation', 0) for c in filtered) / len(filtered), 1) if filtered else 0,
            "total_files": sum(len(c.get('files', [])) for c in filtered)
        }
    })

@app.route('/api/contributions', methods=['POST'])
@rate_limit(max_per_minute=30)
def create_contribution():
    """Enhanced contribution submission"""
    try:
        data = request.json
        
        required_fields = ['molecule_id', 'title', 'type', 'researcher_id']
        for field in required_fields:
            if not data.get(field):
                return jsonify({"success": False, "error": f"Missing field: {field}"}), 400
        
        # Verify molecule exists
        molecule = next((m for m in molecules_db if m["id"] == data['molecule_id']), None)
        if not molecule:
            return jsonify({"success": False, "error": "Molecule not found"}), 404
        
        # Verify researcher exists
        researcher = next((r for r in researchers_db if r["id"] == data['researcher_id']), None)
        if not researcher:
            return jsonify({"success": False, "error": "Researcher not found"}), 404
        
        contribution_id = f"cont_{len(contributions_db) + 1:03d}"
        
        # AI validation of contribution
        ai_validation = random.randint(75, 98)
        reputation_reward = ai_validation // 4
        
        new_contribution = {
            "id": contribution_id,
            "molecule_id": data['molecule_id'],
            "type": data['type'],
            "title": data['title'],
            "researcher_id": data['researcher_id'],
            "description": data.get('description', ''),
            "methodology": data.get('methodology', ''),
            "results": data.get('results', ''),
            "data_url": data.get('data_url', ''),
            "files": data.get('files', []),
            "ai_validation": ai_validation,
            "timestamp": datetime.now().isoformat() + "Z",
            "blockchain_tx": f"0x{hashlib.sha256(contribution_id.encode()).hexdigest()[:40]}...",
            "reputation_reward": reputation_reward,
            "status": "pending_review"
        }
        
        contributions_db.append(new_contribution)
        
        # Update researcher reputation
        researcher['reputation_score'] = min(100, researcher.get('reputation_score', 0) + reputation_reward)
        researcher['total_contributions'] = researcher.get('total_contributions', 0) + 1
        
        # Update molecule contribution count
        molecule['contribution_count'] = molecule.get('contribution_count', 0) + 1
        
        # Add to knowledge graph
        knowledge_graph.nodes.append({
            'id': contribution_id,
            'type': 'contribution',
            'label': data['title'],
            'properties': {'type': data['type']}
        })
        knowledge_graph.edges.append({
            'from': data['researcher_id'],
            'to': contribution_id,
            'relationship': 'CREATED',
            'weight': 1.0
        })
        knowledge_graph.edges.append({
            'from': contribution_id,
            'to': data['molecule_id'],
            'relationship': 'ABOUT',
            'weight': 1.0
        })
        
        return jsonify({
            "success": True,
            "message": "Contribution submitted successfully",
            "contribution": new_contribution,
            "reputation_update": {
                "researcher": researcher['name'],
                "new_score": researcher['reputation_score'],
                "reward": reputation_reward
            },
            "next_steps": [
                "AI validation in progress",
                "Community review available in 24h",
                "Blockchain confirmation pending"
            ]
        }), 201
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/knowledge-graph', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_knowledge_graph():
    """Get knowledge graph data"""
    query = request.args.get('query', 'all')
    
    if query == 'all':
        return jsonify({
            "success": True,
            "graph": {
                "nodes": knowledge_graph.nodes,
                "edges": knowledge_graph.edges
            },
            "statistics": {
                "total_nodes": len(knowledge_graph.nodes),
                "total_edges": len(knowledge_graph.edges),
                "node_types": list(set(n['type'] for n in knowledge_graph.nodes)),
                "relationship_types": list(set(e['relationship'] for e in knowledge_graph.edges))
            }
        })
    
    elif query == 'molecules':
        molecule_nodes = [n for n in knowledge_graph.nodes if n['type'] == 'molecule']
        molecule_edges = [e for e in knowledge_graph.edges 
                         if any(n['id'] == e['from'] or n['id'] == e['to'] 
                               for n in molecule_nodes)]
        
        return jsonify({
            "success": True,
            "subgraph": {
                "nodes": molecule_nodes,
                "edges": molecule_edges
            }
        })
    
    elif query == 'researchers':
        researcher_nodes = [n for n in knowledge_graph.nodes if n['type'] == 'researcher']
        researcher_edges = [e for e in knowledge_graph.edges 
                          if any(n['id'] == e['from'] or n['id'] == e['to'] 
                                for n in researcher_nodes)]
        
        return jsonify({
            "success": True,
            "subgraph": {
                "nodes": researcher_nodes,
                "edges": researcher_edges
            }
        })
    
    return jsonify({"success": False, "error": "Invalid query"}), 400
    

# Additional helper endpoints backed by `database` wrapper
@app.route('/api/search', methods=['GET'])
def search():
    query = request.args.get('q', '')
    category = request.args.get('category', 'all')
    
    results = database.search(query, category)
    return jsonify({
        "success": True,
        "query": query,
        "results": results
    })


@app.route('/api/statistics', methods=['GET'])
def get_statistics():
    return jsonify({
        "success": True,
        "statistics": database.get_statistics()
    })

@app.route('/api/molecules/<molecule_id>/ecosystem', methods=['GET'])
def get_molecule_ecosystem(molecule_id):
    ecosystem = database.get_molecule_ecosystem(molecule_id)
    if ecosystem:
        return jsonify({
            "success": True,
            "ecosystem": ecosystem
        })
    return jsonify({"success": False, "error": "Molecule not found"}), 404
    

@app.route('/api/knowledge-graph/query', methods=['POST'])
@rate_limit(max_per_minute=30)
def query_knowledge_graph():
    """Query knowledge graph with specific parameters"""
    try:
        data = request.json
        query_type = data.get('type')
        
        if query_type == 'find_similar':
            molecule_id = data.get('molecule_id')
            similar = knowledge_graph.query('find_similar', molecule_id=molecule_id)
            
            # Get details of similar molecules
            similar_details = []
            for sim in similar:
                molecule = next((m for m in molecules_db if m['id'] == sim['molecule_id']), None)
                if molecule:
                    similar_details.append({
                        "molecule": molecule,
                        "similarity": sim['similarity']
                    })
            
            return jsonify({
                "success": True,
                "query": f"Find similar to {molecule_id}",
                "results": similar_details,
                "count": len(similar_details)
            })
        
        elif query_type == 'researcher_network':
            researcher_id = data.get('researcher_id')
            network = knowledge_graph.query('researcher_network', researcher_id=researcher_id)
            
            # Get researcher details
            researcher = next((r for r in researchers_db if r['id'] == researcher_id), None)
            
            return jsonify({
                "success": True,
                "query": f"Network for {researcher_id}",
                "researcher": researcher,
                "contributions": network,
                "contribution_count": len(network)
            })
        
        elif query_type == 'molecule_properties':
            molecule_id = data.get('molecule_id')
            properties = knowledge_graph.query('molecule_properties', molecule_id=molecule_id)
            
            return jsonify({
                "success": True,
                "query": f"Properties of {molecule_id}",
                "properties": properties
            })
        
        return jsonify({"success": False, "error": "Invalid query type"}), 400
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/web/search', methods=['POST'])
@rate_limit(max_per_minute=20)
def web_search():
    """Enhanced web scraping/search endpoint"""
    try:
        data = request.json
        query = data.get('query')
        sources = data.get('sources', ['pubchem', 'pubmed', 'clinical_trials'])
        
        if not query:
            return jsonify({"success": False, "error": "Missing search query"}), 400
        
        print(f"🌐 Web agent searching for: {query}")
        
        # Run web search
        job = web_agent.search_all(query)
        
        # Auto-validate against existing molecules
        validation_results = []
        for molecule in molecules_db:
            if query.lower() in molecule['name'].lower():
                validation_results.append({
                    "molecule_id": molecule['id'],
                    "name": molecule['name'],
                    "match_type": "name_match",
                    "confidence": 0.9
                })
            elif query in molecule.get('smiles', ''):
                validation_results.append({
                    "molecule_id": molecule['id'],
                    "name": molecule['name'],
                    "match_type": "smiles_match",
                    "confidence": 1.0
                })
        
        return jsonify({
            "success": True,
            "job_id": job['id'],
            "query": query,
            "sources_searched": list(job['results'].keys()),
            "results": job['results'],
            "auto_validation": {
                "molecules_found": validation_results,
                "total_matches": len(validation_results)
            },
            "timestamp": job['timestamp'],
            "next_actions": [
                "Compare with existing database",
                "Validate with AI agents",
                "Add to knowledge graph"
            ]
        })
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/web/jobs', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_web_jobs():
    """Get web scraping job history"""
    sort_by = request.args.get('sort', 'recent')
    limit = int(request.args.get('limit', 20))
    
    sorted_jobs = web_scraping_jobs.copy()
    if sort_by == 'recent':
        sorted_jobs.sort(key=lambda x: x.get('timestamp', ''), reverse=True)
    
    return jsonify({
        "success": True,
        "jobs": sorted_jobs[:limit],
        "total_jobs": len(web_scraping_jobs),
        "stats": {
            "total_queries": len(web_scraping_jobs),
            "success_rate": 0.95,
            "avg_response_time": "2.3s",
            "cache_hit_rate": 0.65
        }
    })

@app.route('/api/researchers', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_researchers():
    """Get all researchers"""
    sort_by = request.args.get('sort', 'reputation')
    
    sorted_researchers = researchers_db.copy()
    if sort_by == 'reputation':
        sorted_researchers.sort(key=lambda x: x.get('reputation_score', 0), reverse=True)
    elif sort_by == 'contributions':
        sorted_researchers.sort(key=lambda x: x.get('total_contributions', 0), reverse=True)
    elif sort_by == 'name':
        sorted_researchers.sort(key=lambda x: x.get('name', ''))
    
    return jsonify({
        "success": True,
        "count": len(sorted_researchers),
        "researchers": sorted_researchers,
        "leaderboard": [
            {
                "rank": i + 1,
                "name": r['name'],
                "reputation": r['reputation_score'],
                "contributions": r['total_contributions'],
                "institution": r['institution']
            }
            for i, r in enumerate(sorted_researchers[:10])
        ]
    })

@app.route('/api/researchers/<researcher_id>', methods=['GET'])
@rate_limit(max_per_minute=60)
def get_researcher(researcher_id):
    """Get detailed researcher information"""
    researcher = next((r for r in researchers_db if r["id"] == researcher_id), None)
    if not researcher:
        return jsonify({"success": False, "error": "Researcher not found"}), 404
    
    # Get researcher's contributions
    researcher_contributions = [c for c in contributions_db if c["researcher_id"] == researcher_id]
    
    # Get researcher's molecules
    researcher_molecules = []
    for contribution in researcher_contributions:
        molecule = next((m for m in molecules_db if m["id"] == contribution["molecule_id"]), None)
        if molecule and molecule not in researcher_molecules:
            researcher_molecules.append(molecule)
    
    # Get knowledge graph network
    network = knowledge_graph.query('researcher_network', researcher_id=researcher_id)
    
    return jsonify({
        "success": True,
        "researcher": researcher,
        "contributions": {
            "total": len(researcher_contributions),
            "list": researcher_contributions[:10],
            "types": list(set(c['type'] for c in researcher_contributions))
        },
        "molecules": {
            "total": len(researcher_molecules),
            "list": [{"id": m['id'], "name": m['name'], "ai_score": m['ai_score']} 
                    for m in researcher_molecules[:10]]
        },
        "network": {
            "connections": network,
            "collaboration_score": round(len(network) * 0.1, 1)
        },
        "analytics": {
            "avg_ai_validation": round(
                sum(c.get('ai_validation', 0) for c in researcher_contributions) / 
                len(researcher_contributions), 1
            ) if researcher_contributions else 0,
            "total_reputation_earned": sum(c.get('reputation_reward', 0) for c in researcher_contributions),
            "active_since": min(c.get('timestamp', '') for c in researcher_contributions) if researcher_contributions else None
        }
    })

@app.route('/api/dashboard', methods=['GET'])
@rate_limit(max_per_minute=30)
def get_dashboard():
    """Comprehensive dashboard endpoint"""
    # Recent activity
    recent_molecules = sorted(
        molecules_db, 
        key=lambda x: x.get('created_at', ''), 
        reverse=True
    )[:5]
    
    recent_contributions = sorted(
        contributions_db,
        key=lambda x: x.get('timestamp', ''),
        reverse=True
    )[:5]
    
    # Top performers
    top_molecules = sorted(
        molecules_db,
        key=lambda x: x.get('ai_score', 0),
        reverse=True
    )[:5]
    
    top_researchers = sorted(
        researchers_db,
        key=lambda x: x.get('reputation_score', 0),
        reverse=True
    )[:5]
    
    # System metrics
    total_ai_score = sum(m.get('ai_score', 0) for m in molecules_db)
    avg_ai_score = round(total_ai_score / len(molecules_db)) if molecules_db else 0
    
    total_value = sum(
        m.get('blockchain_data', {}).get('value_usd', 0) 
        for m in molecules_db 
        if m.get('on_blockchain')
    )
    
    return jsonify({
        "success": True,
        "dashboard": {
            "overview": {
                "total_molecules": len(molecules_db),
                "verified_molecules": len([m for m in molecules_db if m.get('on_blockchain')]),
                "total_researchers": len(researchers_db),
                "total_contributions": len(contributions_db),
                "total_value_locked": round(total_value, 2),
                "avg_ai_score": avg_ai_score
            },
            "recent_activity": {
                "molecules": recent_molecules,
                "contributions": recent_contributions
            },
            "top_performers": {
                "molecules": top_molecules,
                "researchers": top_researchers
            },
            "system_health": {
                "ai_agents": f"{len([a for a in ai_agents_db if a['status'] == 'active'])}/{len(ai_agents_db)} active",
                "response_time": "0.45s avg",
                "uptime": "99.9%",
                "cache_hit_rate": "92%"
            },
            "blockchain": {
                "flare_price": flare_simulator.ftso_price,
                "total_transactions": len(blockchain_txs),
                "network": "Coston Testnet",
                "status": "Connected"
            }
        },
        "timestamp": datetime.now().isoformat() + "Z",
        "refresh_interval": 30
    })

@app.route('/api/export/molecules', methods=['GET'])
def export_molecules():
    """Export molecules in various formats"""
    format_type = request.args.get('format', 'json')
    
    if format_type == 'json':
        return jsonify({
            "success": True,
            "format": "json",
            "count": len(molecules_db),
            "data": molecules_db,
            "timestamp": datetime.now().isoformat() + "Z"
        })
    
    elif format_type == 'csv':
        # Simple CSV conversion
        import csv
        from io import StringIO
        
        output = StringIO()
        writer = csv.writer(output)
        
        # Write header
        if molecules_db:
            headers = molecules_db[0].keys()
            writer.writerow(headers)
            
            # Write data
            for molecule in molecules_db:
                writer.writerow([str(molecule.get(h, '')) for h in headers])
        
        return jsonify({
            "success": True,
            "format": "csv",
            "data": output.getvalue(),
            "timestamp": datetime.now().isoformat() + "Z"
        })
    
    return jsonify({"success": False, "error": "Unsupported format"}), 400

@app.route('/api/molecules/<molecule_id>/structure', methods=['GET'])
def get_molecule_structure(molecule_id):
    """Get 3D structure data for molecule"""
    molecule = next((m for m in molecules_db if m["id"] == molecule_id), None)
    if not molecule:
        return jsonify({"success": False, "error": "Molecule not found"}), 404
    
    # Generate mock 3D coordinates based on SMILES
    smiles = molecule.get('smiles', '')
    
    # Simple 3D structure generation (for demo)
    structure = {
        "molecule_id": molecule_id,
        "name": molecule.get('name'),
        "format": "xyz",
        "atom_count": len(smiles.replace('=', '').replace('#', '')) // 2,  # Rough estimate
        "coordinates": generate_mock_bonds(smiles),
        "bonds": generate_mock_bonds(smiles),
        "properties": {
            "unit_cell": [10.0, 10.0, 10.0, 90, 90, 90],
            "space_group": "P1",
            "symmetry": "Triclinic"
        }
    }
    
    return jsonify({
        "success": True,
        "structure": structure,
        "visualization_url": f"/api/visualize/{molecule_id}/3d",
        "download_formats": ["xyz", "pdb", "cif", "mol2"]
    })
    


def generate_mock_bonds(smiles):
    """Generate mock bonds for demo"""
    bonds = []
    atom_count = len([c for c in smiles if c.isalpha()])

    for i in range(max(0, atom_count - 1)):
        bonds.append({
            "from": i,
            "to": i + 1,
            "order": 1,
            "length": round(random.uniform(1.0, 1.5), 3)
        })

    return bonds

# ===== ERROR HANDLERS =====
@app.errorhandler(404)
def not_found(error):
    return jsonify({
        "success": False,
        "error": "Endpoint not found",
        "message": "The requested API endpoint does not exist",
        "documentation": "/"
    }), 404

@app.errorhandler(405)
def method_not_allowed(error):
    return jsonify({
        "success": False,
        "error": "Method not allowed",
        "message": "The HTTP method is not supported for this endpoint"
    }), 405

@app.errorhandler(500)
def internal_error(error):
    return jsonify({
        "success": False,
        "error": "Internal server error",
        "message": "An unexpected error occurred"
    }), 500

# ===== ROOT ENDPOINT WITH DOCUMENTATION =====
@app.route('/')
def home():
    return """
    <html>
    <head>
        <title>🧬 MolChain Backend v2.0</title>
        <style>
            body { 
                background: #0f172a; 
                color: white; 
                font-family: 'Segoe UI', system-ui, sans-serif; 
                padding: 40px;
                max-width: 1200px;
                margin: 0 auto;
                line-height: 1.6;
            }
            h1 { 
                color: #60a5fa; 
                font-size: 2.5rem;
                margin-bottom: 10px;
                background: linear-gradient(90deg, #60a5fa, #8b5cf6);
                -webkit-background-clip: text;
                background-clip: text;
                color: transparent;
            }
            .subtitle {
                font-size: 1.2rem;
                opacity: 0.9;
                margin-bottom: 30px;
            }
            .badge {
                display: inline-block;
                padding: 4px 12px;
                background: linear-gradient(90deg, #3b82f6, #8b5cf6);
                color: white;
                border-radius: 20px;
                font-size: 0.9rem;
                font-weight: 600;
                margin: 0 5px 10px 0;
            }
            .endpoint { 
                background: #1e293b; 
                padding: 20px; 
                margin: 15px 0; 
                border-radius: 12px;
                border-left: 4px solid #3b82f6;
                transition: transform 0.2s;
            }
            .endpoint:hover {
                transform: translateY(-2px);
                box-shadow: 0 10px 20px rgba(0, 0, 0, 0.2);
            }
            .method { 
                display: inline-block; 
                padding: 6px 12px; 
                background: #3b82f6; 
                color: white; 
                border-radius: 6px;
                font-weight: bold;
                margin-right: 12px;
                font-family: monospace;
                min-width: 60px;
                text-align: center;
            }
            .method.get { background: #10b981; }
            .method.post { background: #f59e0b; }
            .method.put { background: #8b5cf6; }
            .method.delete { background: #ef4444; }
            .url { 
                color: #cbd5e1; 
                font-family: 'Courier New', monospace;
                font-size: 1.1rem;
            }
            .description {
                margin: 10px 0 0 72px;
                color: #94a3b8;
            }
            .section {
                margin: 40px 0 20px 0;
                padding-bottom: 10px;
                border-bottom: 2px solid #334155;
            }
            .section h2 {
                color: #cbd5e1;
                font-size: 1.5rem;
            }
            .stats {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin: 30px 0;
            }
            .stat-card {
                background: #1e293b;
                padding: 20px;
                border-radius: 10px;
                text-align: center;
                border: 1px solid #334155;
            }
            .stat-value {
                font-size: 2rem;
                font-weight: 800;
                color: #60a5fa;
                margin-bottom: 5px;
            }
            .stat-label {
                font-size: 0.9rem;
                color: #94a3b8;
                text-transform: uppercase;
                letter-spacing: 1px;
            }
            .feature-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                gap: 20px;
                margin: 30px 0;
            }
            .feature-card {
                background: #1e293b;
                padding: 20px;
                border-radius: 10px;
                border: 1px solid #334155;
            }
            .feature-card h3 {
                color: #cbd5e1;
                margin-top: 0;
            }
            .code-block {
                background: #0f172a;
                padding: 15px;
                border-radius: 8px;
                border: 1px solid #334155;
                font-family: 'Courier New', monospace;
                margin: 20px 0;
                overflow-x: auto;
            }
            .footer {
                margin-top: 40px;
                padding-top: 20px;
                border-top: 1px solid #334155;
                text-align: center;
                color: #94a3b8;
                font-size: 0.9rem;
            }
        </style>
    </head>
    <body>
        <h1>🧬 MolChain AI + Blockchain Backend</h1>
        <div class="subtitle">
            Complete backend for ETH Oxford 2026 Hackathon • AI Agent Swarms + Flare Blockchain Integration
        </div>
        
        <div>
            <span class="badge">Flare Integration</span>
            <span class="badge">AI Agent Swarm</span>
            <span class="badge">Knowledge Graph</span>
            <span class="badge">Web Scraping</span>
            <span class="badge">REST API</span>
            <span class="badge">Real-time</span>
        </div>
        
        <div class="stats">
            <div class="stat-card">
                <div class="stat-value">{molecules}</div>
                <div class="stat-label">Molecules</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{agents}</div>
                <div class="stat-label">AI Agents</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{contributions}</div>
                <div class="stat-label">Contributions</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">${value}</div>
                <div class="stat-label">Total Value</div>
            </div>
        </div>
        
        <div class="feature-grid">
            <div class="feature-card">
                <h3>🤖 AI Agent Swarm</h3>
                <p>5 specialized AI agents analyzing molecules concurrently with consensus mechanism.</p>
                <ul>
                    <li>Toxicity Predictor</li>
                    <li>Solubility Analyzer</li>
                    <li>Binding Affinity AI</li>
                    <li>Drug Likeness Evaluator</li>
                    <li>Synthesis Feasibility</li>
                </ul>
            </div>
            
            <div class="feature-card">
                <h3>⛓️ Flare Blockchain</h3>
                <p>Complete Flare protocol integration including FTSO, FDC, and smart contracts.</p>
                <ul>
                    <li>FTSO for molecular valuation</li>
                    <li>Smart contract registry</li>
                    <li>Real-time price feeds</li>
                    <li>On-chain validation</li>
                </ul>
            </div>
            
            <div class="feature-card">
                <h3>🧠 Knowledge Graph</h3>
                <p>Intelligent graph database for molecular relationships and research discovery.</p>
                <ul>
                    <li>Molecular similarity</li>
                    <li>Research networks</li>
                    <li>Property relationships</li>
                    <li>Graph queries</li>
                </ul>
            </div>
        </div>
        
        <div class="section">
            <h2>🚀 Quick Start</h2>
            <div class="code-block">
# Test the backend
curl http://localhost:5000/api/health

# Get all molecules
curl http://localhost:5000/api/molecules

# Run AI analysis
curl -X POST http://localhost:5000/api/analyze \\
  -H "Content-Type: application/json" \\
  -d '{"molecule_id": "mol_001"}'

# Get Flare price
curl http://localhost:5000/api/flare/price
            </div>
        </div>
        
        <div class="section">
            <h2>📡 Core API Endpoints</h2>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/health</span>
                <div class="description">Health check with system status and statistics</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/molecules</span>
                <div class="description">Get all molecules with filtering and pagination</div>
            </div>
            
            <div class="endpoint">
                <span class="method post">POST</span>
                <span class="url">/api/analyze</span>
                <div class="description">Run AI swarm analysis on a molecule</div>
            </div>
            
            <div class="endpoint">
                <span class="method post">POST</span>
                <span class="url">/api/analyze/batch</span>
                <div class="description">Batch analysis of multiple molecules</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/stats</span>
                <div class="description">Comprehensive platform statistics</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/agents</span>
                <div class="description">AI agent swarm status and performance</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/flare/price</span>
                <div class="description">Flare FTSO price feed simulation</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/flare/history</span>
                <div class="description">Flare price history data</div>
            </div>
            
            <div class="endpoint">
                <span class="method post">POST</span>
                <span class="url">/api/flare/register</span>
                <div class="description">Register molecule on Flare blockchain</div>
            </div>
            
            <div class="endpoint">
                <span class="method post">POST</span>
                <span class="url">/api/web/search</span>
                <div class="description">Web scraping agent for molecular data</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/knowledge-graph</span>
                <div class="description">Knowledge graph data and queries</div>
            </div>
            
            <div class="endpoint">
                <span class="method get">GET</span>
                <span class="url">/api/dashboard</span>
                <div class="description">Comprehensive dashboard data</div>
            </div>
        </div>
        
        <div class="section">
            <h2>🔗 React Frontend Integration</h2>
            <div class="code-block">
// React fetch example
const fetchMolecules = async () => {
  const response = await fetch('http://localhost:5000/api/molecules');
  const data = await response.json();
  return data.molecules;
};

// Run AI analysis
const analyzeMolecule = async (moleculeId) => {
  const response = await fetch('http://localhost:5000/api/analyze', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({ molecule_id: moleculeId })
  });
  return await response.json();
};

// Get real-time stats
const getStats = async () => {
  const response = await fetch('http://localhost:5000/api/stats');
  return await response.json();
};
            </div>
        </div>
        
        <div class="footer">
            <p>🧬 MolChain • ETH Oxford 2026 Hackathon Submission</p>
            <p>Backend Server running on port 5000 • Frontend should connect to http://localhost:5000</p>
            <p>Press Ctrl+C to stop the server</p>
        </div>
        
        <script>
            // Update stats with real data
            fetch('/api/health')
                .then(r => r.json())
                .then(data => {
                    const stats = data.stats;
                    document.querySelectorAll('.stat-value')[0].textContent = stats.molecules;
                    document.querySelectorAll('.stat-value')[1].textContent = stats.ai_agents;
                    document.querySelectorAll('.stat-value')[2].textContent = stats.contributions;
                    document.querySelectorAll('.stat-value')[3].textContent = '$' + 
                        Math.round(stats.blockchain_txs * 124.5).toLocaleString();
                });
        </script>
    </body>
    </html>
    """.format(
        molecules=len(molecules_db),
        agents=len(ai_agents_db),
        contributions=len(contributions_db),
        value=round(len(blockchain_txs) * 124.5)
    )

# ===== START SERVER =====
if __name__ == '__main__':
    print("\n" + "="*70)
    print("🚀 MOLCHAIN AI + BLOCKCHAIN BACKEND v2.0 STARTING...")
    print("="*70)
    print("\n🤖 AI Agent Swarm Initialized:")
    for agent_id, agent in ai_swarm.agents.items():
        print(f"  • {agent['name']}: ✓ Ready (Accuracy: {agent['accuracy']*100}%)")
    
    print("\n⛓️ Flare Blockchain Simulator:")
    print(f"  • Network: Coston Testnet (Chain ID: {flare_simulator.chain_id})")
    print(f"  • FTSO Price: ${flare_simulator.ftso_price}")
    print(f"  • Gas Price: {flare_simulator.gas_price / 1e9} Gwei")
    
    print("\n🧠 Knowledge Graph:")
    print(f"  • Nodes: {len(knowledge_graph.nodes)}")
    print(f"  • Edges: {len(knowledge_graph.edges)}")
    
    print("\n📡 Essential Endpoints:")
    print("  Health Check:      http://localhost:5000/api/health")
    print("  Dashboard:         http://localhost:5000/api/dashboard")
    print("  All Molecules:     http://localhost:5000/api/molecules")
    print("  AI Analysis:       POST http://localhost:5000/api/analyze")
    print("  Flare Price:       http://localhost:5000/api/flare/price")
    print("  Knowledge Graph:   http://localhost:5000/api/knowledge-graph")
    print("  Web Search:        POST http://localhost:5000/api/web/search")
    
    print("\n" + "="*70)
    print("✅ Server ready! React frontend should connect automatically")
    print("🎮 Demo Interface: http://localhost:5173")
    print("⏰ Press Ctrl+C to stop the server")
    print("="*70 + "\n")
    
    # Run the server
    app.run(
        debug=True, 
        port=5000, 
        use_reloader=False,
        threaded=True
    )