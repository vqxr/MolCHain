#!/usr/bin/env python3
"""
Simple test script to verify the AI analysis and report download endpoints work
"""
import requests
import json
from datetime import datetime

# Backend URL
BASE_URL = "http://localhost:5000"

def test_create_molecule():
    """Test creating a molecule"""
    print("\n1️⃣  Creating a test molecule...")
    payload = {
        "name": "Aspirin-Test",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "formula": "C9H8O4",
        "molecular_weight": 180.16,
        "description": "Test aspirin molecule for analysis"
    }
    
    response = requests.post(f"{BASE_URL}/api/molecules", json=payload)
    if response.status_code == 201:
        data = response.json()
        if data.get('success'):
            molecule_id = data['molecule']['id']
            print(f"   ✓ Molecule created: {molecule_id}")
            return molecule_id
    
    print(f"   ✗ Failed to create molecule: {response.text}")
    return None

def test_analyze_molecule(molecule_id):
    """Test AI analysis"""
    print(f"\n2️⃣  Running AI analysis on molecule {molecule_id}...")
    payload = {"molecule_id": molecule_id}
    
    response = requests.post(f"{BASE_URL}/api/analyze", json=payload)
    if response.status_code == 200:
        data = response.json()
        if data.get('success'):
            print(f"   ✓ Analysis complete!")
            print(f"   - Consensus Score: {data['consensus_score']}")
            print(f"   - Blockchain Ready: {data['blockchain_ready']}")
            print(f"   - Agents analyzed:")
            for agent_id, result in data['agent_results'].items():
                print(f"     • {result['agent_name']}: {result['score']} (confidence: {result['confidence']})")
            return True
    
    print(f"   ✗ Analysis failed: {response.text}")
    return False

def test_download_report(molecule_id):
    """Test report download"""
    print(f"\n3️⃣  Downloading analysis report for {molecule_id}...")
    payload = {"molecule_id": molecule_id}
    
    response = requests.post(f"{BASE_URL}/api/analyze/report", json=payload)
    if response.status_code == 200:
        # Save the report
        filename = f"test_report_{molecule_id}.txt"
        with open(filename, 'w') as f:
            f.write(response.text)
        print(f"   ✓ Report downloaded and saved as {filename}")
        print(f"   - Report size: {len(response.text)} bytes")
        print(f"   - First 200 chars:")
        print(f"     {response.text[:200]}...")
        return True
    
    print(f"   ✗ Report download failed: {response.status_code}")
    return False

def main():
    print("=" * 70)
    print("   MolChain AI Analysis & Report Download Test")
    print("=" * 70)
    
    # Test 1: Create molecule
    molecule_id = test_create_molecule()
    if not molecule_id:
        print("\n❌ Test failed at molecule creation")
        return
    
    # Test 2: Analyze molecule
    if not test_analyze_molecule(molecule_id):
        print("\n❌ Test failed at analysis")
        return
    
    # Test 3: Download report
    if not test_download_report(molecule_id):
        print("\n❌ Test failed at report download")
        return
    
    print("\n" + "=" * 70)
    print("✅ All tests passed!")
    print("=" * 70)

if __name__ == "__main__":
    main()
