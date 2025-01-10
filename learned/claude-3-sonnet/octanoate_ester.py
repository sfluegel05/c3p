"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:52247 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the acid component is octanoic acid (caprylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group pattern
    ester_pattern = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[OX2][#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # More specific pattern for octanoate group:
    # - Exactly 8 carbons in chain (7 + carbonyl carbon)
    # - No branching allowed on the chain
    # - Must end in ester group
    octanoate_pattern = Chem.MolFromSmarts("[CH3X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CX3](=[OX1])[OX2]")
    
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No linear octanoate chain found"
    
    matches = mol.GetSubstructMatches(octanoate_pattern)
    
    for match in matches:
        # Check that none of the carbons in the chain are part of a ring
        is_in_ring = False
        for atom_idx in match[:-2]:  # Exclude the ester oxygen and carbonyl oxygen
            if mol.GetAtomWithIdx(atom_idx).IsInRing():
                is_in_ring = True
                break
        
        if is_in_ring:
            continue
            
        # Check that the carbons in the chain don't have additional connections
        has_extra_connections = False
        for atom_idx in match[:-2]:  # Exclude the ester group
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetDegree() > 2 and atom_idx != 0:  # First carbon (CH3) can only have 1 connection
                has_extra_connections = True
                break
                
        if has_extra_connections:
            continue
            
        return True, "Contains octanoate ester group (linear 8-carbon chain with terminal ester)"
    
    return False, "No valid octanoate ester group found"

def test_examples():
    """Test function with some known examples"""
    test_cases = [
        ("CCCCCCCC(=O)OCC", True),  # ethyl octanoate
        ("CCCCCCCC(=O)OC", True),   # methyl octanoate
        ("CCCCCCC(=O)OCC", False),  # ethyl heptanoate (too short)
        ("CCCCCCCCC(=O)OCC", False),# ethyl nonanoate (too long)
        ("CCCCCCCCOC", False),      # octyl ether (not an ester)
        ("O(CCCCCCCC)C(=O)CCCCCCC", False),  # octyl octanoate (wrong orientation)
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_octanoate_ester(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()