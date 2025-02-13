"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde has an aldehyde group (-CHO) directly attached to an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_arenecarbaldehyde, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all aldehyde groups
    # Match both explicit and implicit hydrogen representations
    aldehyde_pattern = Chem.MolFromSmarts("[$([CH](=O))]")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not matches:
        return False, "No aldehyde group found"

    # For each aldehyde carbon, check if it's connected to an aromatic system
    for match in matches:
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])
        
        # Get all neighboring atoms of the aldehyde carbon
        for neighbor in aldehyde_carbon.GetNeighbors():
            # Skip the oxygen of the aldehyde group
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:
                continue
                
            # Check if the neighbor is part of an aromatic system
            if neighbor.GetIsAromatic():
                return True, "Contains aldehyde group directly attached to aromatic ring"
            
            # Check if neighbor is part of a ring that might be aromatic
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
            
            for ring in atom_rings:
                if neighbor.GetIdx() in ring:
                    # Check if any atom in this ring is aromatic
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    if any(atom.GetIsAromatic() for atom in ring_atoms):
                        return True, "Contains aldehyde group directly attached to aromatic ring system"

    return False, "No aldehyde group attached to aromatic ring"

def test_smiles():
    """Test function with known examples"""
    test_cases = [
        ("O=CC1=CC=C(N)C=C1", True),  # p-aminobenzaldehyde
        ("Cc1ccccc1C=O", True),        # o-tolualdehyde
        ("O=CC=1N(CC)C=CC1", True),    # 1-Ethyl-1H-pyrrole-2-carboxaldehyde
        ("CCN(CC)c1ccc(C=O)cc1", True),# 4-(diethylamino)benzaldehyde
        ("CC(=O)CC(=O)CC", False),     # Not an arenecarbaldehyde
        ("c1ccccc1", False),           # Not an aldehyde
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_arenecarbaldehyde(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_smiles()