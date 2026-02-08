from rdkit import Chem
from rdkit.Chem import Descriptors

def validate_smiles(smiles):
    """Check if SMILES is valid"""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def simple_analysis(smiles):
    """Simple molecular analysis for hackathon"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}
    
    return {
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Descriptors.MolLogP(mol), 2),
        "h_bond_donors": Descriptors.NumHDonors(mol),
        "h_bond_acceptors": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "ring_count": mol.GetRingInfo().NumRings()
    }