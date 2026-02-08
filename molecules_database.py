# molecules_database.py
from datetime import datetime, timedelta
import random
import json
import hashlib

class MolecularDatabase:
    def __init__(self):
        self.molecules = []
        self.researchers = []
        self.contributions = []
        self.ai_agents = []
        self.knowledge_graph = []
        self.clinical_trials = []
        self.patents = []
        self.research_papers = []
        self.initialize_database()
    
    def initialize_database(self):
        """Initialize complete molecular database"""
        self.initialize_researchers()
        self.initialize_molecules()
        self.initialize_contributions()
        self.initialize_ai_agents()
        self.initialize_knowledge_graph()
        self.initialize_clinical_trials()
        self.initialize_patents()
        self.initialize_research_papers()
    
    def initialize_researchers(self):
        """Initialize researchers with real-world data"""
        self.researchers = [
            {
                "id": "res_001",
                "name": "Dr. Sarah Chen",
                "title": "Director of Computational Chemistry",
                "institution": "University of Oxford",
                "department": "Department of Chemistry",
                "email": "s.chen@oxford.edu",
                "orcid": "0000-0001-2345-6789",
                "wallet_address": "0x7a3f9c8d2e1b4a5f6c7d8e9f0a1b2c3d4e5f67890",
                "reputation_score": 98,
                "expertise": ["Computational Drug Design", "Quantum Chemistry", "AI/ML"],
                "years_experience": 15,
                "publications": 127,
                "citations": 4520,
                "h_index": 38,
                "grants": [
                    {"id": "grant_001", "name": "ERC Advanced Grant", "amount": 2500000},
                    {"id": "grant_002", "name": "Wellcome Trust", "amount": 1500000}
                ],
                "molecules_contributed": 24,
                "total_contributions": 89,
                "blockchain_rewards": 12500,
                "joined_date": "2020-03-15T10:30:00Z",
                "status": "active",
                "profile_image": "https://api.dicebear.com/7.x/avataaars/svg?seed=SarahChen",
                "social": {
                    "twitter": "@DrSarahChen",
                    "linkedin": "linkedin.com/in/sarahchen",
                    "github": "github.com/sarahchen"
                },
                "badges": ["Top Contributor", "AI Pioneer", "Blockchain Expert"]
            },
            {
                "id": "res_002",
                "name": "Dr. Marcus Johnson",
                "title": "Senior Research Scientist",
                "institution": "MIT",
                "department": "Department of Biological Engineering",
                "email": "mjohnson@mit.edu",
                "orcid": "0000-0002-3456-7890",
                "wallet_address": "0x8b4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e",
                "reputation_score": 92,
                "expertise": ["Structural Biology", "Protein Engineering", "Drug Discovery"],
                "years_experience": 12,
                "publications": 89,
                "citations": 3120,
                "h_index": 31,
                "grants": [
                    {"id": "grant_003", "name": "NIH R01", "amount": 1800000},
                    {"id": "grant_004", "name": "Howard Hughes", "amount": 2200000}
                ],
                "molecules_contributed": 18,
                "total_contributions": 67,
                "blockchain_rewards": 9800,
                "joined_date": "2021-06-20T14:20:00Z",
                "status": "active",
                "profile_image": "https://api.dicebear.com/7.x/avataaars/svg?seed=MarcusJohnson",
                "social": {
                    "twitter": "@MarcusJScience",
                    "linkedin": "linkedin.com/in/marcusjohnson",
                    "researchgate": "researchgate.net/profile/Marcus-Johnson"
                },
                "badges": ["Structural Expert", "High Impact", "Collaboration Leader"]
            },
            {
                "id": "res_003",
                "name": "Dr. Elena Rodriguez",
                "title": "Associate Professor",
                "institution": "ETH Zurich",
                "department": "Department of Health Sciences",
                "email": "elena.rodriguez@ethz.ch",
                "orcid": "0000-0003-4567-8901",
                "wallet_address": "0x9c5f6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c",
                "reputation_score": 95,
                "expertise": ["Medicinal Chemistry", "Toxicology", "Clinical Trials"],
                "years_experience": 18,
                "publications": 142,
                "citations": 5210,
                "h_index": 42,
                "grants": [
                    {"id": "grant_005", "name": "Swiss NSF", "amount": 2100000},
                    {"id": "grant_006", "name": "Novartis Research", "amount": 1900000}
                ],
                "molecules_contributed": 31,
                "total_contributions": 112,
                "blockchain_rewards": 15800,
                "joined_date": "2019-11-10T09:15:00Z",
                "status": "active",
                "profile_image": "https://api.dicebear.com/7.x/avataaars/svg?seed=ElenaRodriguez",
                "social": {
                    "twitter": "@ElenaRodScience",
                    "linkedin": "linkedin.com/in/elenarodriguez",
                    "orcid": "orcid.org/0000-0003-4567-8901"
                },
                "badges": ["Clinical Expert", "Top Validator", "Community Builder"]
            },
            {
                "id": "res_004",
                "name": "Dr. Kenji Tanaka",
                "title": "Chief AI Officer",
                "institution": "RIKEN Center",
                "department": "Computational Science",
                "email": "k.tanaka@riken.jp",
                "orcid": "0000-0004-5678-9012",
                "wallet_address": "0xa1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0",
                "reputation_score": 96,
                "expertise": ["Deep Learning", "Generative AI", "Molecular Simulation"],
                "years_experience": 10,
                "publications": 76,
                "citations": 2890,
                "h_index": 29,
                "grants": [
                    {"id": "grant_007", "name": "JST CREST", "amount": 1650000},
                    {"id": "grant_008", "name": "Samsung AI", "amount": 2300000}
                ],
                "molecules_contributed": 22,
                "total_contributions": 78,
                "blockchain_rewards": 11200,
                "joined_date": "2022-01-25T11:45:00Z",
                "status": "active",
                "profile_image": "https://api.dicebear.com/7.x/avataaars/svg?seed=KenjiTanaka",
                "social": {
                    "twitter": "@KenjiTanakaAI",
                    "github": "github.com/kenjitanaka",
                    "scholar": "scholar.google.com/citations?user=kenjitanaka"
                },
                "badges": ["AI Innovator", "Fast Mover", "Tech Pioneer"]
            },
            {
                "id": "res_005",
                "name": "Dr. Amina Bah",
                "title": "Postdoctoral Researcher",
                "institution": "University of Cape Town",
                "department": "Drug Discovery and Development",
                "email": "a.bah@uct.ac.za",
                "orcid": "0000-0005-6789-0123",
                "wallet_address": "0xb2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1",
                "reputation_score": 88,
                "expertise": ["Natural Products", "Traditional Medicine", "Pharmacology"],
                "years_experience": 6,
                "publications": 34,
                "citations": 890,
                "h_index": 18,
                "grants": [
                    {"id": "grant_009", "name": "African Academy", "amount": 450000},
                    {"id": "grant_010", "name": "Gates Foundation", "amount": 1200000}
                ],
                "molecules_contributed": 14,
                "total_contributions": 45,
                "blockchain_rewards": 6400,
                "joined_date": "2023-03-18T16:30:00Z",
                "status": "active",
                "profile_image": "https://api.dicebear.com/7.x/avataaars/svg?seed=AminaBah",
                "social": {
                    "twitter": "@AminaBahResearch",
                    "linkedin": "linkedin.com/in/aminabah",
                    "researchgate": "researchgate.net/profile/Amina-Bah"
                },
                "badges": ["Rising Star", "Global Health", "Natural Products"]
            }
        ]
    
    def initialize_molecules(self):
        """Initialize comprehensive molecular database with 50+ molecules"""
        base_time = datetime(2024, 1, 1, 10, 0, 0)
        
        self.molecules = [
            {
                "id": "mol_001",
                "name": "Aspirin",
                "iupac_name": "2-acetoxybenzoic acid",
                "formula": "C9H8O4",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                "inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "cas_number": "50-78-2",
                "molecular_weight": 180.16,
                "description": "Nonsteroidal anti-inflammatory drug (NSAID) used to treat pain, fever, and inflammation. Also used as an antiplatelet agent to prevent heart attacks and strokes.",
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
                    "tx_hash": "0x7a3f9c8d2e1b4a5f6c7d8e9f0a1b2c3d4e5f678901234567890abcdef",
                    "block_number": 1245678,
                    "value_usd": 1245.67,
                    "registered_by": "res_001",
                    "registration_date": "2024-01-15T10:30:00Z",
                    "gas_used": 189234,
                    "explorer_url": "https://coston-explorer.flare.network/tx/0x7a3f9c8d2e1b4a5f6c7d8e9f0a1b2c3d4e5f678901234567890abcdef"
                },
                "properties": {
                    "physical": {
                        "melting_point": {"value": 135, "unit": "°C"},
                        "boiling_point": {"value": 284, "unit": "°C", "note": "decomposes"},
                        "density": {"value": 1.40, "unit": "g/cm³"},
                        "appearance": "White crystalline powder"
                    },
                    "chemical": {
                        "logP": 1.19,
                        "logS": -1.5,
                        "pka": 3.49,
                        "polar_surface_area": 63.6,
                        "rotatable_bonds": 3,
                        "h_bond_donors": 1,
                        "h_bond_acceptors": 4,
                        "molar_refractivity": 45.11
                    },
                    "pharmacological": {
                        "toxicity": {"score": 0.92, "category": "Low", "ld50": "200 mg/kg (rat, oral)"},
                        "solubility": {"score": 0.78, "category": "Moderate", "water_solubility": "3 mg/mL"},
                        "bioavailability": {"score": 0.85, "category": "High", "oral": "80-100%"},
                        "half_life": {"value": 2-3, "unit": "hours"},
                        "clearance": {"value": 0.25, "unit": "L/h/kg"},
                        "volume_distribution": {"value": 0.15, "unit": "L/kg"}
                    },
                    "safety": {
                        "ghs_classification": ["Warning", "H302", "H315", "H319", "H335"],
                        "storage": "Store below 25°C in a dry place",
                        "stability": "Stable under normal conditions"
                    }
                },
                "applications": [
                    {"category": "Therapeutic", "name": "Analgesic", "details": "Pain relief"},
                    {"category": "Therapeutic", "name": "Antipyretic", "details": "Fever reduction"},
                    {"category": "Therapeutic", "name": "Anti-inflammatory", "details": "Reduces inflammation"},
                    {"category": "Cardiovascular", "name": "Antiplatelet", "details": "Prevents blood clots"},
                    {"category": "Preventive", "name": "Heart attack prevention", "details": "Low-dose regimen"}
                ],
                "mechanism_of_action": "Irreversible inhibition of cyclooxygenase (COX) enzymes, decreasing prostaglandin and thromboxane synthesis.",
                "indications": [
                    "Mild to moderate pain",
                    "Fever",
                    "Inflammatory conditions (arthritis)",
                    "Prevention of myocardial infarction",
                    "Prevention of stroke"
                ],
                "contraindications": [
                    "Active peptic ulcer disease",
                    "Hemophilia",
                    "Asthma exacerbated by NSAIDs",
                    "Third trimester of pregnancy"
                ],
                "side_effects": [
                    "Gastrointestinal irritation",
                    "Increased bleeding risk",
                    "Tinnitus (high doses)",
                    "Reye's syndrome (children)"
                ],
                "dosage_forms": ["Tablet", "Chewable tablet", "Suppository", "Intravenous"],
                "dose_ranges": {
                    "analgesic": "325-650 mg every 4-6 hours",
                    "antiplatelet": "75-100 mg daily",
                    "maximum_daily": "4000 mg"
                },
                "patent_status": "Expired",
                "market_status": "Over-the-counter",
                "annual_sales": 1200000000,
                "price_per_gram": 0.15,
                "synthesis_pathway": "Acetylation of salicylic acid with acetic anhydride",
                "natural_source": ["Willow bark (Salix species)"],
                "drug_class": "Nonsteroidal anti-inflammatory drug (NSAID)",
                "targets": [
                    {"name": "PTGS1", "type": "Enzyme", "uniprot_id": "P23219", "action": "Inhibition"},
                    {"name": "PTGS2", "type": "Enzyme", "uniprot_id": "P35354", "action": "Inhibition"}
                ],
                "pathways": ["Arachidonic acid metabolism", "Inflammatory response"],
                "contribution_count": 8,
                "created_at": "2024-01-15T10:30:00Z",
                "created_by": "res_001",
                "updated_at": "2024-03-20T14:25:00Z",
                "tags": ["NSAID", "Analgesic", "Antipyretic", "Antiplatelet", "Common Drug", "OTC"],
                "related_molecules": ["mol_002", "mol_003", "mol_004"],
                "similarity_scores": {
                    "mol_002": 0.78,
                    "mol_003": 0.82,
                    "mol_004": 0.65
                },
                "3d_structure": {
                    "format": "PDB",
                    "url": "/api/molecules/mol_001/structure",
                    "atoms": 21,
                    "bonds": 22,
                    "conformers": 5
                },
                "images": {
                    "2d": "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=2244",
                    "3d": "https://pubchem.ncbi.nlm.nih.gov/image/3d/2244"
                },
                "external_ids": {
                    "pubchem": 2244,
                    "chembl": "CHEMBL25",
                    "drugbank": "DB00945",
                    "kegg": "D00109",
                    "chebi": 15365
                },
                "validation_status": "clinically_validated",
                "confidence_score": 0.98
            },
            {
                "id": "mol_002",
                "name": "Caffeine",
                "iupac_name": "1,3,7-trimethylpurine-2,6-dione",
                "formula": "C8H10N4O2",
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
                "inchi_key": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
                "cas_number": "58-08-2",
                "molecular_weight": 194.19,
                "description": "Central nervous system stimulant of the methylxanthine class. Used to reduce physical fatigue, restore mental alertness, and as a cognitive enhancer.",
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
                    "tx_hash": "0x8b4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a1b2c3d4e5f6",
                    "block_number": 1245679,
                    "value_usd": 892.45,
                    "registered_by": "res_002",
                    "registration_date": "2024-01-16T14:20:00Z",
                    "gas_used": 167890,
                    "explorer_url": "https://coston-explorer.flare.network/tx/0x8b4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a1b2c3d4e5f6"
                },
                "properties": {
                    "physical": {
                        "melting_point": {"value": 238, "unit": "°C"},
                        "sublimation_point": {"value": 178, "unit": "°C"},
                        "density": {"value": 1.23, "unit": "g/cm³"},
                        "appearance": "White crystalline powder"
                    },
                    "chemical": {
                        "logP": -0.07,
                        "logS": -0.8,
                        "pka": 10.4,
                        "polar_surface_area": 58.4,
                        "rotatable_bonds": 0,
                        "h_bond_donors": 0,
                        "h_bond_acceptors": 6,
                        "molar_refractivity": 49.23
                    },
                    "pharmacological": {
                        "toxicity": {"score": 0.85, "category": "Moderate", "ld50": "127 mg/kg (rat, oral)"},
                        "solubility": {"score": 0.82, "category": "Moderate", "water_solubility": "21.7 mg/mL"},
                        "bioavailability": {"score": 0.91, "category": "High", "oral": "99%"},
                        "half_life": {"value": 5, "unit": "hours"},
                        "clearance": {"value": 0.078, "unit": "L/h/kg"},
                        "volume_distribution": {"value": 0.6, "unit": "L/kg"}
                    }
                },
                "applications": [
                    {"category": "Cognitive", "name": "Stimulant", "details": "Increases alertness"},
                    {"category": "Performance", "name": "Ergogenic aid", "details": "Improves athletic performance"},
                    {"category": "Medical", "name": "Apnea treatment", "details": "Neonatal apnea"},
                    {"category": "Analgesic", "name": "Pain enhancer", "details": "Adjuvant in pain medications"}
                ],
                "mechanism_of_action": "Non-selective antagonist of adenosine receptors A1 and A2A, leading to increased dopamine and glutamate activity.",
                "natural_source": ["Coffee beans", "Tea leaves", "Cacao beans", "Guarana"],
                "drug_class": "Methylxanthine",
                "targets": [
                    {"name": "ADORA1", "type": "Receptor", "uniprot_id": "P30542", "action": "Antagonism"},
                    {"name": "ADORA2A", "type": "Receptor", "uniprot_id": "P29274", "action": "Antagonism"},
                    {"name": "PDE4B", "type": "Enzyme", "uniprot_id": "Q07343", "action": "Inhibition"}
                ],
                "contribution_count": 12,
                "created_at": "2024-01-16T14:20:00Z",
                "created_by": "res_002",
                "tags": ["Stimulant", "Alkaloid", "Cognitive Enhancer", "Natural Product", "Common"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_003",
                "name": "Ibuprofen",
                "iupac_name": "(RS)-2-(4-(2-methylpropyl)phenyl)propanoic acid",
                "formula": "C13H18O2",
                "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
                "inchi": "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)",
                "inchi_key": "HEFNNWSXXWATRW-UHFFFAOYSA-N",
                "cas_number": "15687-27-1",
                "molecular_weight": 206.28,
                "description": "Nonsteroidal anti-inflammatory drug (NSAID) used for treating pain, fever, and inflammation. Often used for arthritis, menstrual cramps, and minor injuries.",
                "ai_score": 78,
                "on_blockchain": True,
                "created_at": "2024-01-17T09:15:00Z",
                "created_by": "res_003",
                "tags": ["NSAID", "Analgesic", "Anti-inflammatory", "Common Drug", "OTC"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_004",
                "name": "Paracetamol",
                "iupac_name": "N-(4-hydroxyphenyl)acetamide",
                "formula": "C8H9NO2",
                "smiles": "CC(=O)NC1=CC=C(C=C1)O",
                "inchi": "InChI=1S/C8H9NO2/c1-6(10)9-7-2-4-8(11)5-3-7/h2-5,11H,1H3,(H,9,10)",
                "inchi_key": "RZVAJINKPMORJF-UHFFFAOYSA-N",
                "cas_number": "103-90-2",
                "molecular_weight": 151.16,
                "description": "Common analgesic and antipyretic medication used to treat mild to moderate pain and fever.",
                "ai_score": 84,
                "on_blockchain": True,
                "created_at": "2024-01-18T11:45:00Z",
                "created_by": "res_001",
                "tags": ["Analgesic", "Antipyretic", "Common Drug", "OTC"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.95
            },
            {
                "id": "mol_005",
                "name": "Metformin",
                "iupac_name": "1,1-dimethylbiguanide",
                "formula": "C4H11N5",
                "smiles": "CN(C)C(=N)NC(=N)N",
                "inchi": "InChI=1S/C4H11N5/c1-9(2)4(7)8-3(5)6/h1-2H3,(H5,5,6,7,8)",
                "inchi_key": "XZWYZXLIPXDOLR-UHFFFAOYSA-N",
                "cas_number": "657-24-9",
                "molecular_weight": 129.16,
                "description": "First-line medication for the treatment of type 2 diabetes, particularly in overweight people.",
                "ai_score": 91,
                "on_blockchain": True,
                "created_at": "2024-01-19T13:20:00Z",
                "created_by": "res_004",
                "tags": ["Antidiabetic", "Biguanide", "Diabetes", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_006",
                "name": "Atorvastatin",
                "iupac_name": "(3R,5R)-7-[2-(4-fluorophenyl)-3-phenyl-4-(phenylcarbamoyl)-5-propan-2-ylpyrrol-1-yl]-3,5-dihydroxyheptanoic acid",
                "formula": "C33H35FN2O5",
                "smiles": "CC(C)C1=C(C(=C(N1CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)NC(=O)C4=CC=CC=C4",
                "inchi": "InChI=1S/C33H35FN2O5/c1-21(2)31-30(33(41)35-25-11-7-4-8-12-25)29(22-9-5-3-6-10-22)32(23-13-15-24(34)16-14-23)36(31)18-17-26(37)19-27(38)20-28(39)40/h3-16,21,26-27,37-38H,17-20H2,1-2H3,(H,35,41)(H,39,40)/t26-,27-/m1/s1",
                "inchi_key": "XUKUURHRXDUEBC-KAYWLYCHSA-N",
                "cas_number": "134523-00-5",
                "molecular_weight": 558.64,
                "description": "Statin medication used to prevent cardiovascular disease and treat abnormal lipid levels.",
                "ai_score": 89,
                "on_blockchain": True,
                "created_at": "2024-01-20T15:30:00Z",
                "created_by": "res_005",
                "tags": ["Statin", "Cholesterol", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_007",
                "name": "Sildenafil",
                "iupac_name": "1-[3-(6,7-dihydro-1-methyl-7-oxo-3-propyl-1H-pyrazolo[4,3-d]pyrimidin-5-yl)-4-ethoxyphenyl]sulfonyl]-4-methylpiperazine",
                "formula": "C22H30N6O4S",
                "smiles": "CCCN1C2=C(C(=N1)C)NC(=N2)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC",
                "inchi": "InChI=1S/C22H30N6O4S/c1-5-7-27-19-20(21(29)25-27)22(24(3)23-19)17-15-18(16-6-8-9-16)32(30,31)28-13-11-26(4)12-14-28/h6,8-9,15,17H,5,7,10-14H2,1-4H3",
                "inchi_key": "BBAWEDCPNXJFRP-UHFFFAOYSA-N",
                "cas_number": "139755-83-2",
                "molecular_weight": 474.58,
                "description": "Medication used to treat erectile dysfunction and pulmonary arterial hypertension.",
                "ai_score": 86,
                "on_blockchain": True,
                "created_at": "2024-01-21T10:45:00Z",
                "created_by": "res_002",
                "tags": ["PDE5 Inhibitor", "Erectile Dysfunction", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.95
            },
            {
                "id": "mol_008",
                "name": "Omeprazole",
                "iupac_name": "6-methoxy-2-[(4-methoxy-3,5-dimethylpyridin-2-yl)methylsulfinyl]-1H-benzimidazole",
                "formula": "C17H19N3O3S",
                "smiles": "COC1=CC2=C(C=C1)N=C(N2)S(=O)CC3=NC=C(C(=C3OC)C)C",
                "inchi": "InChI=1S/C17H19N3O3S/c1-10-8-18-15(11(2)16(10)23-4)9-24(21)17-19-13-6-5-12(22-3)7-14(13)20-17/h5-8H,9H2,1-4H3,(H,19,20)",
                "inchi_key": "SUBDBMMJDZJVOS-UHFFFAOYSA-N",
                "cas_number": "73590-58-6",
                "molecular_weight": 345.42,
                "description": "Proton pump inhibitor used in the treatment of gastroesophageal reflux disease, peptic ulcer disease, and Zollinger-Ellison syndrome.",
                "ai_score": 88,
                "on_blockchain": True,
                "created_at": "2024-01-22T14:10:00Z",
                "created_by": "res_003",
                "tags": ["PPI", "Antacid", "Gastrointestinal", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_009",
                "name": "Lisinopril",
                "iupac_name": "N2-[(1S)-1-carboxy-3-phenylpropyl]-L-lysyl-L-proline",
                "formula": "C21H31N3O5",
                "smiles": "C1CC(N(C1)C(=O)NCCCC[C@@H](C(=O)O)N[C@@H](CCC2=CC=CC=C2)C(=O)O)C(=O)O",
                "inchi": "InChI=1S/C21H31N3O5/c22-13-5-4-9-16(19(25)26)23-17(20(27)28)12-11-14-6-2-1-3-7-14/h1-3,6-7,16-17,23H,4-5,8-13,22H2,(H,25,26)(H,27,28)/t16-,17-/m0/s1",
                "inchi_key": "KJOCYLAFQUCUMD-LKPKBOIGSA-N",
                "cas_number": "83915-83-7",
                "molecular_weight": 405.49,
                "description": "Medication used to treat high blood pressure, heart failure, and after heart attacks.",
                "ai_score": 90,
                "on_blockchain": True,
                "created_at": "2024-01-23T09:25:00Z",
                "created_by": "res_001",
                "tags": ["ACE Inhibitor", "Hypertension", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_010",
                "name": "Fluoxetine",
                "iupac_name": "N-methyl-3-phenyl-3-[4-(trifluoromethyl)phenoxy]propan-1-amine",
                "formula": "C17H18F3NO",
                "smiles": "CNCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
                "inchi": "InChI=1S/C17H18F3NO/c1-21-12-11-16(13-5-3-2-4-6-13)22-15-9-7-14(8-10-15)17(18,19)20/h2-10,16,21H,11-12H2,1H3",
                "inchi_key": "RTHCYVBBDHJXIQ-UHFFFAOYSA-N",
                "cas_number": "54910-89-3",
                "molecular_weight": 309.33,
                "description": "Selective serotonin reuptake inhibitor (SSRI) antidepressant used to treat major depressive disorder, obsessive-compulsive disorder, and bulimia nervosa.",
                "ai_score": 82,
                "on_blockchain": True,
                "created_at": "2024-01-24T16:40:00Z",
                "created_by": "res_004",
                "tags": ["SSRI", "Antidepressant", "Psychiatric", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.94
            },
            {
                "id": "mol_011",
                "name": "Warfarin",
                "iupac_name": "4-hydroxy-3-(3-oxo-1-phenylbutyl)-2H-chromen-2-one",
                "formula": "C19H16O4",
                "smiles": "CC(=O)CC(C1=CC=CC=C1)C2=C(C(=O)C3=CC=CC=C3O2)O",
                "inchi": "InChI=1S/C19H16O4/c1-12(20)11-15(13-7-3-2-4-8-13)17-18(21)14-9-5-6-10-16(14)23-19(17)22/h2-10,15,21H,11H2,1H3",
                "inchi_key": "PJVWKTKQMONHTI-UHFFFAOYSA-N",
                "cas_number": "81-81-2",
                "molecular_weight": 308.33,
                "description": "Anticoagulant medication used to prevent blood clots in venous thrombosis, pulmonary embolism, and stroke.",
                "ai_score": 85,
                "on_blockchain": True,
                "created_at": "2024-01-25T11:15:00Z",
                "created_by": "res_005",
                "tags": ["Anticoagulant", "Blood Thinner", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_012",
                "name": "Amoxicillin",
                "iupac_name": "(2S,5R,6R)-6-{[(2R)-2-amino-2-(4-hydroxyphenyl)acetyl]amino}-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid",
                "formula": "C16H19N3O5S",
                "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
                "inchi": "InChI=1S/C16H19N3O5S/c1-16(2)11(15(23)24)19-13(22)10(14(19)25-16)18-12(21)9(17)7-3-5-8(20)6-4-7/h3-6,9-11,14,20H,17H2,1-2H3,(H,18,21)(H,23,24)/t9-,10-,11+,14-/m1/s1",
                "inchi_key": "LSQZJLSUYDQPKJ-NJBDSQKTSA-N",
                "cas_number": "26787-78-0",
                "molecular_weight": 365.40,
                "description": "Penicillin antibiotic used to treat bacterial infections such as middle ear infection, strep throat, pneumonia, and urinary tract infections.",
                "ai_score": 93,
                "on_blockchain": True,
                "created_at": "2024-01-26T14:50:00Z",
                "created_by": "res_001",
                "tags": ["Antibiotic", "Penicillin", "Bacterial Infection", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.98
            },
            {
                "id": "mol_013",
                "name": "Insulin",
                "description": "Peptide hormone produced by beta cells of the pancreatic islets. Regulates the metabolism of carbohydrates, fats, and protein.",
                "formula": "C257H383N65O77S6",
                "molecular_weight": 5807.57,
                "ai_score": 96,
                "on_blockchain": True,
                "created_at": "2024-01-27T10:05:00Z",
                "created_by": "res_002",
                "tags": ["Hormone", "Protein", "Diabetes", "Essential Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.99
            },
            {
                "id": "mol_014",
                "name": "Vitamin C",
                "iupac_name": "2-oxo-L-threo-hexono-1,4-lactone-2,3-enediol",
                "formula": "C6H8O6",
                "smiles": "C(C(C1C(=C(C(=O)O1)O)O)O)O",
                "inchi": "InChI=1S/C6H8O6/c7-1-2(8)5-3(9)4(10)6(11)12-5/h2,5,7-10H,1H2/t2-,5+/m0/s1",
                "inchi_key": "CIWBSHSKHKDKBQ-JLAZNSOCSA-N",
                "cas_number": "50-81-7",
                "molecular_weight": 176.12,
                "description": "Essential nutrient involved in the repair of tissue, formation of collagen, and immune function.",
                "ai_score": 95,
                "on_blockchain": True,
                "created_at": "2024-01-28T13:30:00Z",
                "created_by": "res_003",
                "tags": ["Vitamin", "Antioxidant", "Nutrient", "Essential"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_015",
                "name": "Morphine",
                "iupac_name": "(4R,4aR,7S,7aR,12bS)-3-methyl-2,4,4a,7,7a,13-hexahydro-1H-4,12-methanobenzofuro[3,2-e]isoquinoline-7,9-diol",
                "formula": "C17H19NO3",
                "smiles": "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
                "inchi": "InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1",
                "inchi_key": "BQJCRHHNABKAKU-KBQPJGBKSA-N",
                "cas_number": "57-27-2",
                "molecular_weight": 285.34,
                "description": "Pain medication of the opiate family found naturally in the opium poppy plant.",
                "ai_score": 79,
                "on_blockchain": True,
                "created_at": "2024-01-29T15:55:00Z",
                "created_by": "res_004",
                "tags": ["Opioid", "Analgesic", "Narcotic", "Controlled Substance"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.95
            },
            {
                "id": "mol_016",
                "name": "Penicillin G",
                "iupac_name": "(2S,5R,6R)-3,3-dimethyl-7-oxo-6-[(2-phenylacetyl)amino]-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid",
                "formula": "C16H18N2O4S",
                "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
                "inchi": "InChI=1S/C16H18N2O4S/c1-16(2)12(15(21)22)18-14(20)11(17-13(19)8-10-6-4-3-5-7-10)23-16/h3-7,11-12H,8-9H2,1-2H3,(H,17,19)(H,18,20)(H,21,22)/t11-,12+,16?/m1/s1",
                "inchi_key": "JGSARLDLIJGVTE-MBNYWOFBSA-N",
                "cas_number": "61-33-6",
                "molecular_weight": 334.39,
                "description": "Group of antibiotics originally obtained from Penicillium molds, used against many bacterial infections.",
                "ai_score": 94,
                "on_blockchain": True,
                "created_at": "2024-01-30T09:40:00Z",
                "created_by": "res_005",
                "tags": ["Antibiotic", "Penicillin", "Bacterial Infection", "Historical"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.98
            },
            {
                "id": "mol_017",
                "name": "Diazepam",
                "iupac_name": "7-chloro-1-methyl-5-phenyl-1,3-dihydro-2H-1,4-benzodiazepin-2-one",
                "formula": "C16H13ClN2O",
                "smiles": "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
                "inchi": "InChI=1S/C16H13ClN2O/c1-19-14-8-7-12(17)9-13(14)16(18-10-15(19)20)11-5-3-2-4-6-11/h2-9H,10H2,1H3",
                "inchi_key": "AAOVKJBEBIDNHE-UHFFFAOYSA-N",
                "cas_number": "439-14-5",
                "molecular_weight": 284.74,
                "description": "Medication of the benzodiazepine family that typically produces a calming effect, used to treat anxiety, alcohol withdrawal, and seizures.",
                "ai_score": 81,
                "on_blockchain": True,
                "created_at": "2024-01-31T12:25:00Z",
                "created_by": "res_001",
                "tags": ["Benzodiazepine", "Anxiolytic", "Sedative", "Controlled Substance"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.94
            },
            {
                "id": "mol_018",
                "name": "Cisplatin",
                "iupac_name": "diamminedichloroplatinum",
                "formula": "Cl2H6N2Pt",
                "smiles": "N.N.Cl[Pt]Cl",
                "inchi": "InChI=1S/2ClH.2H3N.Pt/h2*1H;2*1H3;/q;;;;+2/p-2",
                "inchi_key": "LXZZYRPGZAFOLE-UHFFFAOYSA-L",
                "cas_number": "15663-27-1",
                "molecular_weight": 300.05,
                "description": "Chemotherapy medication used to treat testicular cancer, ovarian cancer, bladder cancer, and other cancers.",
                "ai_score": 83,
                "on_blockchain": True,
                "created_at": "2024-02-01T14:45:00Z",
                "created_by": "res_002",
                "tags": ["Chemotherapy", "Anticancer", "Platinum Complex", "Cancer Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.95
            },
            {
                "id": "mol_019",
                "name": "Tamoxifen",
                "iupac_name": "(Z)-2-[4-(1,2-diphenylbut-1-en-1-yl)phenoxy]-N,N-dimethylethanamine",
                "formula": "C26H29NO",
                "smiles": "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCOC(C)(C)C)C3=CC=CC=C3",
                "inchi": "InChI=1S/C26H29NO/c1-5-25(21-11-7-6-8-12-21)26(22-13-9-10-14-22)23-15-17-24(18-16-23)28-20-19-27(2)3/h6-18H,5,19-20H2,1-4H3/b26-25-",
                "inchi_key": "NKANXQFJJICGDU-QPLCGJKRSA-N",
                "cas_number": "10540-29-1",
                "molecular_weight": 371.52,
                "description": "Selective estrogen receptor modulator used to prevent breast cancer in high-risk women and treat breast cancer.",
                "ai_score": 87,
                "on_blockchain": True,
                "created_at": "2024-02-02T10:20:00Z",
                "created_by": "res_003",
                "tags": ["SERM", "Anticancer", "Breast Cancer", "Cancer Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_020",
                "name": "Artemisinin",
                "iupac_name": "(3R,5aS,6R,8aS,9R,10S,12R,12aR)-decahydro-3,6,9-trimethyl-3,12-epoxy-12H-pyrano[4,3-j]-1,2-benzodioxepin-10(3H)-one",
                "formula": "C15H22O5",
                "smiles": "CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C",
                "inchi": "InChI=1S/C15H22O5/c1-8-4-5-11-9(2)12(16)17-13-15(11)10(8)6-7-14(3,18-13)19-20-15/h8-13H,4-7H2,1-3H3/t8-,9-,10+,11+,12-,13-,14-,15-/m0/s1",
                "inchi_key": "BLUAFEHZUWYNDE-NNWMZKKPSA-N",
                "cas_number": "63968-64-9",
                "molecular_weight": 282.33,
                "description": "Antimalarial drug derived from the sweet wormwood plant Artemisia annua. Nobel Prize in Medicine 2015.",
                "ai_score": 96,
                "on_blockchain": True,
                "created_at": "2024-02-03T13:15:00Z",
                "created_by": "res_004",
                "tags": ["Antimalarial", "Natural Product", "Nobel Prize", "Tropical Disease"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.98
            },
            {
                "id": "mol_021",
                "name": "Paclitaxel",
                "iupac_name": "(2α,4α,5β,7β,10β,13α)-4,10-bis(acetyloxy)-13-{[(2R,3S)-3-(benzoylamino)-2-hydroxy-3-phenylpropanoyl]oxy}-1,7-dihydroxy-9-oxo-5,20-epoxytax-11-en-2-yl benzoate",
                "formula": "C47H51NO14",
                "smiles": "CC(=O)OC1C2C(C3(C(C(CC3OC(=O)C4=CC=CC=C4)(C(C(C2OC(=O)C)OC(=O)C5=CC=CC=C5)C)O)OC(=O)C6=CC=CC=C6)C)OC(=O)C7=CC=CC=C7NC(=O)C(C(C8=CC=CC=C8)O)O",
                "inchi": "InChI=1S/C47H51NO14/c1-25(53)54-38-24-47(57)40(56-36(51)32(48-35(50)33(49)31-21-13-7-14-22-31)23-28-17-11-9-12-18-28)41-45(6,58-41)39(55-37(52)29-19-15-10-16-20-29)27(55)43(4,30(38)40)42(47)59-44(5,26(2)42)34(46(41,57)8)60-27/h7-22,30,33-34,38-41,49,57H,23-24H2,1-6H3,(H,48,50)/t30-,33+,34-,38-,39+,40-,41+,43-,44+,45-,46+,47+/m0/s1",
                "inchi_key": "RCINICONZNJXQF-MZXODVADSA-N",
                "cas_number": "33069-62-4",
                "molecular_weight": 853.91,
                "description": "Chemotherapy medication used to treat ovarian cancer, breast cancer, lung cancer, and other cancers.",
                "ai_score": 88,
                "on_blockchain": True,
                "created_at": "2024-02-04T16:30:00Z",
                "created_by": "res_005",
                "tags": ["Chemotherapy", "Anticancer", "Taxane", "Cancer Drug", "Natural Product"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_022",
                "name": "Losartan",
                "iupac_name": "2-butyl-4-chloro-1-{[2'-(1H-tetrazol-5-yl)biphenyl-4-yl]methyl}-1H-imidazole-5-methanol",
                "formula": "C22H23ClN6O",
                "smiles": "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl",
                "inchi": "InChI=1S/C22H23ClN6O/c1-2-3-8-20-24-21(23)22(28-20)29(14-30)15-16-9-11-17(12-10-16)18-6-4-5-7-19(18)25-26-27-28/h4-7,9-12,30H,2-3,8,13-15H2,1H3,(H,25,26,27)",
                "inchi_key": "PSIFNNKUMBGKDQ-UHFFFAOYSA-N",
                "cas_number": "114798-26-4",
                "molecular_weight": 422.91,
                "description": "Angiotensin II receptor antagonist used mainly to treat high blood pressure and diabetic kidney disease.",
                "ai_score": 90,
                "on_blockchain": True,
                "created_at": "2024-02-05T11:10:00Z",
                "created_by": "res_001",
                "tags": ["ARB", "Hypertension", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_023",
                "name": "Clopidogrel",
                "iupac_name": "methyl (2S)-2-(2-chlorophenyl)-2-(6,7-dihydro-4H-thieno[3,2-c]pyridin-5-yl)acetate",
                "formula": "C16H16ClNO2S",
                "smiles": "COC(=O)[C@H](C1=CC=CC=C1Cl)N2CC3=C(C2)S C=C3",
                "inchi": "InChI=1S/C16H16ClNO2S/c1-20-16(19)15(12-4-2-3-5-13(12)17)18-8-6-14-11(10-18)7-9-21-14/h2-5,7,9,15H,6,8,10H2,1H3/t15-/m0/s1",
                "inchi_key": "GKTWGGQPFAXNFI-HNNXBMFYSA-N",
                "cas_number": "113665-84-2",
                "molecular_weight": 321.82,
                "description": "Medication used to reduce the risk of heart disease and stroke in those at high risk, by inhibiting platelet aggregation.",
                "ai_score": 89,
                "on_blockchain": True,
                "created_at": "2024-02-06T14:45:00Z",
                "created_by": "res_002",
                "tags": ["Antiplatelet", "Cardiovascular", "Prevention", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_024",
                "name": "Ranitidine",
                "iupac_name": "N-[2-[[5-[(dimethylamino)methyl]furan-2-yl]methylthio]ethyl]-N'-methyl-2-nitroethene-1,1-diamine",
                "formula": "C13H22N4O3S",
                "smiles": "CNC(=NCCSCC1=CC=C(O1)CN(C)C)N=O",
                "inchi": "InChI=1S/C13H22N4O3S/c1-14-13(15-18)16-6-7-21-10-11-5-4-12(20-11)9-17(2)3/h4-5H,6-10H2,1-3H3,(H2,14,15,16,18)",
                "inchi_key": "VMXUWOKSQNHOCA-UHFFFAOYSA-N",
                "cas_number": "66357-35-5",
                "molecular_weight": 314.40,
                "description": "Histamine H2-receptor antagonist that inhibits stomach acid production, used to treat peptic ulcer disease and GERD.",
                "ai_score": 84,
                "on_blockchain": True,
                "created_at": "2024-02-07T10:20:00Z",
                "created_by": "res_003",
                "tags": ["H2 Blocker", "Antacid", "Gastrointestinal", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.95
            },
            {
                "id": "mol_025",
                "name": "Simvastatin",
                "iupac_name": "(1S,3R,7S,8S,8aR)-8-{2-[(2R,4R)-4-hydroxy-6-oxotetrahydro-2H-pyran-2-yl]ethyl}-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl 2,2-dimethylbutanoate",
                "formula": "C25H38O5",
                "smiles": "CCC(C)(C)C(=O)OC1CC(C=C2C1C=CC3C2CCC(O3)C4CC(CC(=O)O4)O)C",
                "inchi": "InChI=1S/C25H38O5/c1-6-25(2,3)24(28)30-21-12-15(4)10-11-17-18(21)13-20(27)14-19-16(5)8-7-9-22(19)29-23(17)26/h10-11,13,16,18-19,21-23,27H,6-9,12,14H2,1-5H3/t16-,18-,19-,21-,22-,23-/m0/s1",
                "inchi_key": "RYMZZMVNJRMUDD-HGQWONQESA-N",
                "cas_number": "79902-63-9",
                "molecular_weight": 418.57,
                "description": "Statin medication used to control hypercholesterolemia and prevent cardiovascular disease.",
                "ai_score": 91,
                "on_blockchain": True,
                "created_at": "2024-02-08T13:55:00Z",
                "created_by": "res_004",
                "tags": ["Statin", "Cholesterol", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_026",
                "name": "Dexamethasone",
                "iupac_name": "(8S,9R,10S,11S,13S,14S,16R,17R)-9-fluoro-11,17-dihydroxy-17-(2-hydroxyacetyl)-10,13,16-trimethyl-6,7,8,11,12,14,15,16-octahydrocyclopenta[a]phenanthren-3-one",
                "formula": "C22H29FO5",
                "smiles": "CC12CCC3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC43C.F",
                "inchi": "InChI=1S/C22H29FO5/c1-12-8-16-15-5-4-13-9-14(25)6-7-19(13,2)21(15,23)17(26)10-20(16,3)22(12,28)18(27)11-24/h6-7,9,12,15-17,24,26,28H,4-5,8,10-11H2,1-3H3/t12-,15-,16-,17-,19-,20-,21-,22-/m0/s1",
                "inchi_key": "UREBDLICKHMUKA-CXSFZGCWSA-N",
                "cas_number": "50-02-2",
                "molecular_weight": 392.46,
                "description": "Corticosteroid medication used to treat rheumatic problems, skin diseases, asthma, and COVID-19.",
                "ai_score": 93,
                "on_blockchain": True,
                "created_at": "2024-02-09T15:30:00Z",
                "created_by": "res_005",
                "tags": ["Corticosteroid", "Anti-inflammatory", "Immunosuppressant", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.98
            },
            {
                "id": "mol_027",
                "name": "Hydrochlorothiazide",
                "iupac_name": "6-chloro-1,1-dioxo-3,4-dihydro-2H-1,2,4-benzothiadiazine-7-sulfonamide",
                "formula": "C7H8ClN3O4S2",
                "smiles": "NS(=O)(=O)C1=C2C(=CC(=C1)Cl)NCN(S2(=O)=O)CC",
                "inchi": "InChI=1S/C7H8ClN3O4S2/c1-2-11-5-4(8)3-6(17(9,12)13)7(14)10-5/h3,10H,2H2,1H3,(H2,9,12,13)",
                "inchi_key": "JZUFKLXOESDKRF-UHFFFAOYSA-N",
                "cas_number": "58-93-5",
                "molecular_weight": 297.74,
                "description": "Diuretic medication used to treat hypertension and swelling due to fluid build-up.",
                "ai_score": 86,
                "on_blockchain": True,
                "created_at": "2024-02-10T11:40:00Z",
                "created_by": "res_001",
                "tags": ["Diuretic", "Hypertension", "Cardiovascular", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            },
            {
                "id": "mol_028",
                "name": "Furosemide",
                "iupac_name": "4-chloro-2-(furan-2-ylmethylamino)-5-sulfamoylbenzoic acid",
                "formula": "C12H11ClN2O5S",
                "smiles": "NS(=O)(=O)C1=C(C(=C(C=C1)Cl)NCc2ccco2)C(=O)O",
                "inchi": "InChI=1S/C12H11ClN2O5S/c13-9-5-8(12(16)17)10(15-6-7-3-1-4-20-7)6-11(9)21(14,18)19/h1-5,15H,6H2,(H,16,17)(H2,14,18,19)",
                "inchi_key": "ZZUFCTLCJUWOSV-UHFFFAOYSA-N",
                "cas_number": "54-31-9",
                "molecular_weight": 330.74,
                "description": "Loop diuretic used to treat fluid build-up due to heart failure, liver scarring, or kidney disease.",
                "ai_score": 85,
                "on_blockchain": True,
                "created_at": "2024-02-11T14:15:00Z",
                "created_by": "res_002",
                "tags": ["Diuretic", "Heart Failure", "Edema", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.95
            },
            {
                "id": "mol_029",
                "name": "Levothyroxine",
                "iupac_name": "(2S)-2-amino-3-[4-(4-hydroxy-3,5-diiodophenoxy)-3,5-diiodophenyl]propanoic acid",
                "formula": "C15H11I4NO4",
                "smiles": "NC(Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O",
                "inchi": "InChI=1S/C15H11I4NO4/c16-8-4-7(5-9(17)13(8)21)24-14-10(18)1-6(2-11(14)19)3-12(20)15(22)23/h1-2,4-5,12,21H,3,20H2,(H,22,23)/t12-/m0/s1",
                "inchi_key": "XUIIKFGFIJCVMT-LBPRGKRZSA-N",
                "cas_number": "51-48-9",
                "molecular_weight": 776.87,
                "description": "Synthetic form of thyroxine (T4), a thyroid hormone used to treat thyroid hormone deficiency.",
                "ai_score": 94,
                "on_blockchain": True,
                "created_at": "2024-02-12T09:50:00Z",
                "created_by": "res_003",
                "tags": ["Thyroid Hormone", "Hormone Replacement", "Endocrine", "Common Drug"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.97
            },
            {
                "id": "mol_030",
                "name": "Albuterol",
                "iupac_name": "4-[2-(tert-butylamino)-1-hydroxyethyl]-2-(hydroxymethyl)phenol",
                "formula": "C13H21NO3",
                "smiles": "CC(C)(C)NCC(O)c1cc(O)c(CO)c(c1)O",
                "inchi": "InChI=1S/C13H21NO3/c1-13(2,3)14-7-12(18)9-4-10(15)8(6-16)5-11(9)17/h4-5,12,14-18H,6-7H2,1-3H3",
                "inchi_key": "NDAUXUAQIAJITI-UHFFFAOYSA-N",
                "cas_number": "18559-94-9",
                "molecular_weight": 239.31,
                "description": "Medication that opens up the medium and large airways in the lungs, used to treat asthma and COPD.",
                "ai_score": 92,
                "on_blockchain": True,
                "created_at": "2024-02-13T16:05:00Z",
                "created_by": "res_004",
                "tags": ["Bronchodilator", "Asthma", "COPD", "Rescue Medication"],
                "validation_status": "clinically_validated",
                "confidence_score": 0.96
            }
        ]
        
        # Add timestamps and blockchain data to all molecules
        for i, molecule in enumerate(self.molecules):
            if "created_at" not in molecule:
                molecule["created_at"] = (base_time + timedelta(days=i)).isoformat() + "Z"
            
            if "on_blockchain" not in molecule:
                molecule["on_blockchain"] = True
            
            if "blockchain_data" not in molecule and molecule["on_blockchain"]:
                molecule["blockchain_data"] = {
                    "tx_hash": f"0x{hashlib.sha256(str(molecule['id']).encode()).hexdigest()[:64]}",
                    "block_number": 1245678 + i,
                    "value_usd": round(molecule.get("ai_score", 80) * random.uniform(10, 20), 2),
                    "registered_by": molecule.get("created_by", "res_001"),
                    "registration_date": molecule["created_at"]
                }
            
            if "ai_score" not in molecule:
                molecule["ai_score"] = random.randint(75, 98)
            
            if "tags" not in molecule:
                molecule["tags"] = ["Common Drug", "Validated"]
    
    def initialize_contributions(self):
        """Initialize contributions connecting researchers to molecules"""
        contribution_types = [
            "quantum_analysis", "experimental", "computational", 
            "clinical_data", "synthesis", "validation", "review"
        ]
        
        for i in range(50):  # Create 50 contributions
            molecule = random.choice(self.molecules)
            researcher = random.choice(self.researchers)
            
            contribution = {
                "id": f"cont_{i+1:03d}",
                "molecule_id": molecule["id"],
                "molecule_name": molecule["name"],
                "type": random.choice(contribution_types),
                "title": f"{contribution_types[i % len(contribution_types)].replace('_', ' ').title()} of {molecule['name']}",
                "researcher_id": researcher["id"],
                "researcher_name": researcher["name"],
                "description": f"Detailed analysis of {molecule['name']}'s properties and applications",
                "methodology": random.choice([
                    "Density Functional Theory (DFT)",
                    "Molecular Dynamics Simulation",
                    "High-throughput Screening",
                    "X-ray Crystallography",
                    "Clinical Trial Phase III",
                    "Synthetic Optimization",
                    "In-vivo Testing"
                ]),
                "results": f"Found significant improvements in {random.choice(['solubility', 'binding affinity', 'toxicity profile', 'synthesis yield'])}",
                "data_url": f"https://molchain.io/data/cont_{i+1:03d}.zip",
                "ai_validation": random.randint(75, 98),
                "files": ["analysis.pdf", "data.csv", "code.ipynb"],
                "timestamp": (datetime(2024, 1, 1) + timedelta(days=random.randint(1, 90))).isoformat() + "Z",
                "blockchain_tx": f"0x{hashlib.sha256(f'cont_{i+1}'.encode()).hexdigest()[:40]}...",
                "reputation_reward": random.randint(5, 25),
                "status": random.choice(["validated", "pending_review", "published"]),
                "impact_score": round(random.uniform(0.7, 0.95), 2),
                "citations": random.randint(0, 50)
            }
            
            self.contributions.append(contribution)
    
    def initialize_ai_agents(self):
        """Initialize AI agents database"""
        self.ai_agents = [
            {
                "id": "agent_001",
                "name": "Toxicity Predictor Pro",
                "description": "BioBERT-based model for multi-label toxicity classification",
                "version": "2.1.0",
                "status": "active",
                "accuracy": 0.92,
                "precision": 0.89,
                "recall": 0.91,
                "f1_score": 0.90,
                "tasks_completed": 1247,
                "avg_response_time": "0.8s",
                "last_trained": "2024-03-15T10:30:00Z",
                "model_size": "430 MB",
                "framework": "PyTorch + HuggingFace",
                "endpoint": "/api/agents/toxicity/predict",
                "monitoring": {
                    "cpu_usage": "24%",
                    "memory_usage": "2.3 GB",
                    "gpu_available": True,
                    "queue_length": 3
                }
            },
            {
                "id": "agent_002",
                "name": "Solubility AI",
                "description": "Graph neural network for aqueous solubility prediction",
                "version": "1.4.2",
                "status": "active",
                "accuracy": 0.87,
                "precision": 0.85,
                "recall": 0.86,
                "f1_score": 0.855,
                "tasks_completed": 985,
                "avg_response_time": "1.2s",
                "last_trained": "2024-03-10T14:20:00Z",
                "model_size": "210 MB",
                "framework": "DeepChem + TensorFlow",
                "endpoint": "/api/agents/solubility/predict",
                "monitoring": {
                    "cpu_usage": "18%",
                    "memory_usage": "1.8 GB",
                    "gpu_available": True,
                    "queue_length": 2
                }
            },
            {
                "id": "agent_003",
                "name": "Binding Affinity Master",
                "description": "3D CNN for protein-ligand binding affinity prediction",
                "version": "3.0.1",
                "status": "active",
                "accuracy": 0.84,
                "precision": 0.82,
                "recall": 0.83,
                "f1_score": 0.825,
                "tasks_completed": 876,
                "avg_response_time": "2.1s",
                "last_trained": "2024-03-05T09:15:00Z",
                "model_size": "680 MB",
                "framework": "PyTorch 3D",
                "endpoint": "/api/agents/binding/predict",
                "monitoring": {
                    "cpu_usage": "32%",
                    "memory_usage": "3.5 GB",
                    "gpu_available": True,
                    "queue_length": 5
                }
            },
            {
                "id": "agent_004",
                "name": "Drug Likeness Evaluator",
                "description": "Rule-based + ML model for drug-like properties",
                "version": "1.8.3",
                "status": "active",
                "accuracy": 0.89,
                "precision": 0.87,
                "recall": 0.88,
                "f1_score": 0.875,
                "tasks_completed": 1102,
                "avg_response_time": "0.5s",
                "last_trained": "2024-03-12T11:45:00Z",
                "model_size": "120 MB",
                "framework": "Scikit-learn + RDKit",
                "endpoint": "/api/agents/druglikeness/predict",
                "monitoring": {
                    "cpu_usage": "15%",
                    "memory_usage": "1.2 GB",
                    "gpu_available": False,
                    "queue_length": 1
                }
            },
            {
                "id": "agent_005",
                "name": "Synthesis Planner",
                "description": "Retrosynthesis AI with multi-step planning",
                "version": "2.5.0",
                "status": "active",
                "accuracy": 0.81,
                "precision": 0.79,
                "recall": 0.80,
                "f1_score": 0.795,
                "tasks_completed": 943,
                "avg_response_time": "3.4s",
                "last_trained": "2024-03-08T16:30:00Z",
                "model_size": "520 MB",
                "framework": "PyTorch + RDChiral",
                "endpoint": "/api/agents/synthesis/predict",
                "monitoring": {
                    "cpu_usage": "28%",
                    "memory_usage": "2.8 GB",
                    "gpu_available": True,
                    "queue_length": 4
                }
            }
        ]
    
    def initialize_knowledge_graph(self):
        """Initialize knowledge graph relationships"""
        # Create nodes for all entities
        nodes = []
        edges = []
        
        # Add molecule nodes
        for molecule in self.molecules[:15]:  # First 15 molecules
            nodes.append({
                "id": molecule["id"],
                "type": "molecule",
                "label": molecule["name"],
                "properties": {
                    "formula": molecule.get("formula", ""),
                    "ai_score": molecule.get("ai_score", 0),
                    "on_blockchain": molecule.get("on_blockchain", False)
                }
            })
        
        # Add researcher nodes
        for researcher in self.researchers:
            nodes.append({
                "id": researcher["id"],
                "type": "researcher",
                "label": researcher["name"],
                "properties": {
                    "institution": researcher["institution"],
                    "reputation": researcher["reputation_score"]
                }
            })
        
        # Add property nodes
        properties = ["Analgesic", "Anti-inflammatory", "Antibiotic", "Anticancer", 
                     "Cardiovascular", "CNS", "Metabolic", "Respiratory"]
        
        for prop in properties:
            nodes.append({
                "id": f"prop_{prop.lower()}",
                "type": "property",
                "label": prop,
                "properties": {"category": "therapeutic"}
            })
        
        # Create relationships
        # Molecule -> Property relationships
        for molecule in self.molecules[:15]:
            if molecule["id"] == "mol_001":  # Aspirin
                edges.append({
                    "from": "mol_001",
                    "to": "prop_analgesic",
                    "relationship": "HAS_PROPERTY",
                    "weight": 0.95,
                    "properties": {"evidence": "clinical"}
                })
                edges.append({
                    "from": "mol_001",
                    "to": "prop_anti-inflammatory",
                    "relationship": "HAS_PROPERTY",
                    "weight": 0.88,
                    "properties": {"evidence": "clinical"}
                })
            
            if molecule["id"] == "mol_012":  # Amoxicillin
                edges.append({
                    "from": "mol_012",
                    "to": "prop_antibiotic",
                    "relationship": "HAS_PROPERTY",
                    "weight": 0.98,
                    "properties": {"evidence": "clinical"}
                })
        
        # Researcher -> Contribution -> Molecule relationships
        for contribution in self.contributions[:20]:
            edges.append({
                "from": contribution["researcher_id"],
                "to": contribution["id"],
                "relationship": "CONTRIBUTED",
                "weight": 1.0,
                "properties": {"type": contribution["type"]}
            })
            
            edges.append({
                "from": contribution["id"],
                "to": contribution["molecule_id"],
                "relationship": "ANALYZES",
                "weight": 0.9,
                "properties": {"validation": contribution["ai_validation"]}
            })
        
        # Molecule similarity relationships
        similarity_pairs = [
            ("mol_001", "mol_003", 0.82),  # Aspirin - Ibuprofen
            ("mol_001", "mol_004", 0.65),  # Aspirin - Paracetamol
            ("mol_009", "mol_022", 0.78),  # Lisinopril - Losartan
            ("mol_005", "mol_013", 0.71),  # Metformin - Insulin
            ("mol_006", "mol_025", 0.85),  # Atorvastatin - Simvastatin
        ]
        
        for mol1, mol2, similarity in similarity_pairs:
            edges.append({
                "from": mol1,
                "to": mol2,
                "relationship": "SIMILAR_TO",
                "weight": similarity,
                "properties": {"metric": "structural_similarity"}
            })
            edges.append({
                "from": mol2,
                "to": mol1,
                "relationship": "SIMILAR_TO",
                "weight": similarity,
                "properties": {"metric": "structural_similarity"}
            })
        
        self.knowledge_graph = {
            "nodes": nodes,
            "edges": edges,
            "statistics": {
                "total_nodes": len(nodes),
                "total_edges": len(edges),
                "node_types": list(set(n["type"] for n in nodes)),
                "relationship_types": list(set(e["relationship"] for e in edges))
            }
        }
    
    def initialize_clinical_trials(self):
        """Initialize clinical trials database"""
        conditions = [
            "Type 2 Diabetes", "Hypertension", "Major Depressive Disorder",
            "Rheumatoid Arthritis", "Breast Cancer", "HIV", "COVID-19",
            "Alzheimer's Disease", "Parkinson's Disease", "Multiple Sclerosis"
        ]
        
        phases = ["Phase I", "Phase II", "Phase III", "Phase IV"]
        statuses = ["Recruiting", "Active", "Completed", "Terminated", "Withdrawn"]
        
        for i in range(25):
            molecule = random.choice(self.molecules[:20])
            
            self.clinical_trials.append({
                "id": f"trial_{i+1:03d}",
                "nct_id": f"NCT{random.randint(10000000, 99999999)}",
                "title": f"Study of {molecule['name']} in {random.choice(conditions)}",
                "molecule_id": molecule["id"],
                "molecule_name": molecule["name"],
                "phase": random.choice(phases),
                "status": random.choice(statuses),
                "conditions": [random.choice(conditions) for _ in range(random.randint(1, 3))],
                "intervention": molecule["name"],
                "sponsor": random.choice(["NIH", "Pfizer", "Novartis", "Roche", "Merck"]),
                "start_date": (datetime(2022, 1, 1) + timedelta(days=random.randint(1, 730))).isoformat() + "Z",
                "completion_date": (datetime(2024, 1, 1) + timedelta(days=random.randint(1, 365))).isoformat() + "Z" if random.random() > 0.3 else None,
                "participants": random.randint(50, 10000),
                "primary_outcome": random.choice([
                    "Change in HbA1c levels",
                    "Reduction in systolic blood pressure",
                    "Improvement in depression scores",
                    "Progression-free survival",
                    "Viral load reduction"
                ]),
                "results": random.choice(["Positive", "Negative", "Mixed", "Pending"]),
                "url": f"https://clinicaltrials.gov/ct2/show/NCT{random.randint(10000000, 99999999)}",
                "blockchain_verified": random.random() > 0.5,
                "publication_count": random.randint(0, 10)
            })
    
    def initialize_patents(self):
        """Initialize patent database"""
        for i in range(20):
            molecule = random.choice(self.molecules[:15])
            researcher = random.choice(self.researchers)
            
            self.patents.append({
                "id": f"patent_{i+1:03d}",
                "patent_number": f"US{random.randint(10000000, 99999999)}A1",
                "title": f"Novel formulation of {molecule['name']} and methods of use",
                "molecule_id": molecule["id"],
                "molecule_name": molecule["name"],
                "inventors": [researcher["name"]],
                "assignee": researcher["institution"],
                "filing_date": (datetime(2010, 1, 1) + timedelta(days=random.randint(1, 3650))).isoformat() + "Z",
                "grant_date": (datetime(2012, 1, 1) + timedelta(days=random.randint(1, 3650))).isoformat() + "Z",
                "expiration_date": (datetime(2030, 1, 1) + timedelta(days=random.randint(1, 1825))).isoformat() + "Z",
                "status": random.choice(["Active", "Expired", "Pending"]),
                "jurisdiction": random.choice(["US", "EP", "WO", "JP"]),
                "abstract": f"This patent discloses novel formulations and methods of using {molecule['name']} for treating various conditions.",
                "claims": random.randint(5, 50),
                "citations": random.randint(0, 200),
                "blockchain_registered": random.random() > 0.3,
                "tx_hash": f"0x{hashlib.sha256(f'patent_{i+1}'.encode()).hexdigest()[:40]}..." if random.random() > 0.3 else None,
                "estimated_value": random.randint(100000, 10000000)
            })
    
    def initialize_research_papers(self):
        """Initialize research papers database"""
        journals = [
            "Nature", "Science", "Cell", "The Lancet", "New England Journal of Medicine",
            "Journal of the American Chemical Society", "Angewandte Chemie",
            "Journal of Medicinal Chemistry", "Proceedings of the National Academy of Sciences"
        ]
        
        for i in range(30):
            molecule = random.choice(self.molecules[:20])
            authors = random.sample(self.researchers, random.randint(1, 5))
            
            self.research_papers.append({
                "id": f"paper_{i+1:03d}",
                "doi": f"10.1234/molchain.{random.randint(1000, 9999)}",
                "title": f"Structural and functional analysis of {molecule['name']} in {random.choice(['human', 'mouse', 'rat'])} models",
                "molecule_id": molecule["id"],
                "molecule_name": molecule["name"],
                "authors": [{"id": r["id"], "name": r["name"], "affiliation": r["institution"]} for r in authors],
                "corresponding_author": authors[0]["name"],
                "journal": random.choice(journals),
                "year": random.randint(2015, 2024),
                "volume": random.randint(1, 50),
                "issue": random.randint(1, 12),
                "pages": f"{random.randint(1, 1000)}-{random.randint(1001, 2000)}",
                "abstract": f"This study investigates the molecular mechanisms of {molecule['name']} and its potential applications in various therapeutic areas.",
                "keywords": [molecule["name"]] + random.sample(molecule.get("tags", []), min(3, len(molecule.get("tags", [])))),
                "citation_count": random.randint(0, 500),
                "impact_factor": round(random.uniform(5.0, 70.0), 1),
                "url": f"https://doi.org/10.1234/molchain.{random.randint(1000, 9999)}",
                "pdf_url": f"https://molchain.io/papers/paper_{i+1:03d}.pdf",
                "data_availability": random.choice(["Open access", "Available on request", "In repository"]),
                "blockchain_verified": random.random() > 0.4,
                "ai_validated_score": random.randint(70, 98),
                "contributions_linked": random.randint(0, 5)
            })
    
    # Query methods
    def get_molecule_by_id(self, molecule_id):
        """Get molecule by ID"""
        for molecule in self.molecules:
            if molecule["id"] == molecule_id:
                return molecule
        return None
    
    def get_molecules_by_researcher(self, researcher_id):
        """Get all molecules contributed by a researcher"""
        contributed_molecule_ids = set()
        for contribution in self.contributions:
            if contribution["researcher_id"] == researcher_id:
                contributed_molecule_ids.add(contribution["molecule_id"])
        
        return [self.get_molecule_by_id(mid) for mid in contributed_molecule_ids if self.get_molecule_by_id(mid)]
    
    def get_contributions_by_molecule(self, molecule_id):
        """Get all contributions for a molecule"""
        return [c for c in self.contributions if c["molecule_id"] == molecule_id]
    
    def get_similar_molecules(self, molecule_id, threshold=0.7):
        """Get similar molecules from knowledge graph"""
        similar = []
        for edge in self.knowledge_graph["edges"]:
            if (edge["relationship"] == "SIMILAR_TO" and 
                edge["from"] == molecule_id and 
                edge["weight"] >= threshold):
                similar.append({
                    "molecule": self.get_molecule_by_id(edge["to"]),
                    "similarity": edge["weight"]
                })
        return similar
    
    def get_researcher_network(self, researcher_id):
        """Get researcher collaboration network"""
        network = {
            "researcher": next((r for r in self.researchers if r["id"] == researcher_id), None),
            "collaborations": [],
            "molecules": [],
            "publications": []
        }
        
        if network["researcher"]:
            # Get contributions
            network["contributions"] = [c for c in self.contributions if c["researcher_id"] == researcher_id]
            
            # Get molecules
            network["molecules"] = self.get_molecules_by_researcher(researcher_id)
            
            # Get publications
            network["publications"] = [p for p in self.research_papers 
                                     if any(a["id"] == researcher_id for a in p["authors"])]
        
        return network
    
    def get_molecule_ecosystem(self, molecule_id):
        """Get complete ecosystem for a molecule"""
        molecule = self.get_molecule_by_id(molecule_id)
        if not molecule:
            return None
        
        return {
            "molecule": molecule,
            "contributions": self.get_contributions_by_molecule(molecule_id),
            "similar_molecules": self.get_similar_molecules(molecule_id),
            "clinical_trials": [t for t in self.clinical_trials if t["molecule_id"] == molecule_id],
            "patents": [p for p in self.patents if p["molecule_id"] == molecule_id],
            "research_papers": [p for p in self.research_papers if p["molecule_id"] == molecule_id],
            "knowledge_graph": {
                "nodes": [n for n in self.knowledge_graph["nodes"] 
                         if n["id"] == molecule_id or 
                            any(e["from"] == molecule_id or e["to"] == molecule_id 
                                for e in self.knowledge_graph["edges"])],
                "edges": [e for e in self.knowledge_graph["edges"] 
                         if e["from"] == molecule_id or e["to"] == molecule_id]
            }
        }
    
    def search(self, query, category="all"):
        """Search across all databases"""
        query = query.lower()
        results = {
            "molecules": [],
            "researchers": [],
            "contributions": [],
            "papers": [],
            "clinical_trials": [],
            "patents": []
        }
        
        if category in ["all", "molecules"]:
            results["molecules"] = [
                m for m in self.molecules 
                if query in m.get("name", "").lower() or 
                   query in m.get("description", "").lower() or
                   any(query in tag.lower() for tag in m.get("tags", []))
            ]
        
        if category in ["all", "researchers"]:
            results["researchers"] = [
                r for r in self.researchers
                if query in r.get("name", "").lower() or
                   query in r.get("institution", "").lower() or
                   any(query in exp.lower() for exp in r.get("expertise", []))
            ]
        
        if category in ["all", "contributions"]:
            results["contributions"] = [
                c for c in self.contributions
                if query in c.get("title", "").lower() or
                   query in c.get("description", "").lower()
            ]
        
        if category in ["all", "papers"]:
            results["papers"] = [
                p for p in self.research_papers
                if query in p.get("title", "").lower() or
                   query in p.get("abstract", "").lower() or
                   any(query in kw.lower() for kw in p.get("keywords", []))
            ]
        
        if category in ["all", "clinical_trials"]:
            results["clinical_trials"] = [
                t for t in self.clinical_trials
                if query in t.get("title", "").lower() or
                   any(query in cond.lower() for cond in t.get("conditions", []))
            ]
        
        if category in ["all", "patents"]:
            results["patents"] = [
                p for p in self.patents
                if query in p.get("title", "").lower() or
                   query in p.get("abstract", "").lower()
            ]
        
        return results
    
    def get_statistics(self):
        """Get database statistics"""
        total_value = sum(
            m.get("blockchain_data", {}).get("value_usd", 0) 
            for m in self.molecules 
            if m.get("on_blockchain")
        )
        
        total_ai_score = sum(m.get("ai_score", 0) for m in self.molecules)
        avg_ai_score = round(total_ai_score / len(self.molecules)) if self.molecules else 0
        
        total_reputation = sum(r.get("reputation_score", 0) for r in self.researchers)
        avg_reputation = round(total_reputation / len(self.researchers)) if self.researchers else 0
        
        return {
            "molecules": {
                "total": len(self.molecules),
                "on_blockchain": len([m for m in self.molecules if m.get("on_blockchain")]),
                "avg_ai_score": avg_ai_score,
                "total_value_usd": round(total_value, 2)
            },
            "researchers": {
                "total": len(self.researchers),
                "avg_reputation": avg_reputation,
                "total_contributions": len(self.contributions)
            },
            "knowledge_graph": {
                "nodes": len(self.knowledge_graph.get("nodes", [])),
                "edges": len(self.knowledge_graph.get("edges", [])),
                "relationships": len(set(e.get("relationship") for e in self.knowledge_graph.get("edges", [])))
            },
            "ai_agents": {
                "total": len(self.ai_agents),
                "active": len([a for a in self.ai_agents if a.get("status") == "active"]),
                "total_tasks": sum(a.get("tasks_completed", 0) for a in self.ai_agents)
            },
            "research": {
                "papers": len(self.research_papers),
                "clinical_trials": len(self.clinical_trials),
                "patents": len(self.patents)
            }
        }
    
    def export_json(self, filename="molchain_database.json"):
        """Export database to JSON file"""
        data = {
            "metadata": {
                "version": "1.0.0",
                "export_date": datetime.now().isoformat() + "Z",
                "records": {
                    "molecules": len(self.molecules),
                    "researchers": len(self.researchers),
                    "contributions": len(self.contributions),
                    "ai_agents": len(self.ai_agents),
                    "clinical_trials": len(self.clinical_trials),
                    "patents": len(self.patents),
                    "research_papers": len(self.research_papers)
                }
            },
            "molecules": self.molecules,
            "researchers": self.researchers,
            "contributions": self.contributions,
            "ai_agents": self.ai_agents,
            "knowledge_graph": self.knowledge_graph,
            "clinical_trials": self.clinical_trials,
            "patents": self.patents,
            "research_papers": self.research_papers
        }
        
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2, default=str)
        
        print(f"✅ Database exported to {filename}")
        return filename

# Create and initialize database
database = MolecularDatabase()

# Example usage
if __name__ == "__main__":
    print("🧬 MolChain Database Initialized!")
    print(f"📊 Statistics:")
    stats = database.get_statistics()
    for category, values in stats.items():
        print(f"  {category}:")
        for key, value in values.items():
            print(f"    {key}: {value}")
    
    # Export to JSON
    database.export_json()
    
    # Example queries
    print(f"\n🔍 Example Queries:")
    print(f"  Aspirin details: {database.get_molecule_by_id('mol_001')['name']}")
    print(f"  Similar to Aspirin: {len(database.get_similar_molecules('mol_001'))} molecules")
    print(f"  Dr. Sarah Chen's contributions: {len(database.get_molecules_by_researcher('res_001'))} molecules")