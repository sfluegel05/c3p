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

    # Find all aldehyde groups directly connected to aromatic atoms
    # [CH,CH2]=O attached to aromatic atom
    aldehyde_pattern = Chem.MolFromSmarts("[a]!@[$([CH]=O)]")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if matches:
        # Verify each match
        for match in matches:
            aromatic_atom = mol.GetAtomWithIdx(match[0])
            aldehyde_carbon = mol.GetAtomWithIdx(match[1])
            
            # Double check the aldehyde carbon
            if aldehyde_carbon.GetTotalNumHs() != 1:
                continue
                
            # Verify the aromatic atom is part of a valid aromatic system
            if aromatic_atom.GetIsAromatic():
                # Get the ring this atom belongs to
                ring_info = mol.GetRingInfo()
                atom_rings = ring_info.AtomRings()
                
                for ring in atom_rings:
                    if aromatic_atom.GetIdx() in ring:
                        # Verify the ring is truly aromatic (all atoms in ring are aromatic)
                        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                        if all(atom.GetIsAromatic() for atom in ring_atoms):
                            return True, "Contains aldehyde group directly attached to aromatic ring"

    # Alternative pattern for formyl groups
    formyl_pattern = Chem.MolFromSmarts("[a]!@C(=O)[H]")
    matches = mol.GetSubstructMatches(formyl_pattern)
    
    if matches:
        # Similar verification as above
        for match in matches:
            aromatic_atom = mol.GetAtomWithIdx(match[0])
            formyl_carbon = mol.GetAtomWithIdx(match[1])
            
            if aromatic_atom.GetIsAromatic():
                ring_info = mol.GetRingInfo()
                atom_rings = ring_info.AtomRings()
                
                for ring in atom_rings:
                    if aromatic_atom.GetIdx() in ring:
                        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                        if all(atom.GetIsAromatic() for atom in ring_atoms):
                            return True, "Contains aldehyde group directly attached to aromatic ring"

    return False, "No aldehyde group directly attached to aromatic ring"

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