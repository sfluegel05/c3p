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

    # Look for ester group with 7-carbon chain attached
    # Note: The carbonyl carbon counts as the 8th carbon
    # [CH3] - start with methyl
    # [CH2][CH2][CH2][CH2][CH2][CH2] - six CH2 groups
    # C(=O)O - ester group
    # [#6,#1;!R] means carbon or hydrogen, not in ring (to prevent matching cyclic structures)
    octanoate_pattern = Chem.MolFromSmarts(
        "[CH3][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[OX2]"
    )
    
    if not mol.HasSubstructMatch(octanoate_pattern):
        # Try alternative pattern that allows for substitutions
        alt_pattern = Chem.MolFromSmarts(
            "[CX4][CX4][CX4][CX4][CX4][CX4][CX4]C(=O)[OX2]"
        )
        if not mol.HasSubstructMatch(alt_pattern):
            return False, "No octanoate ester group found"
        
        # Verify the chain isn't part of something longer
        matches = mol.GetSubstructMatches(alt_pattern)
        for match in matches:
            chain_atoms = list(match[:-2])  # Exclude ester oxygen and carbonyl carbon
            
            # Check that first carbon isn't connected to another carbon chain
            first_c = mol.GetAtomWithIdx(chain_atoms[0])
            non_match_neighbors = [n for n in first_c.GetNeighbors() 
                                 if n.GetIdx() not in chain_atoms 
                                 and n.GetAtomicNum() == 6]
            
            if not non_match_neighbors:
                return True, "Contains octanoate ester group with possible substitutions"
                
        return False, "Found similar structure but not a true octanoate ester"

    return True, "Contains octanoate ester group"

def test_examples():
    """Test function with some known examples"""
    test_cases = [
        ("CCCCCCCC(=O)OCC", True),  # ethyl octanoate
        ("CCCCCCCC(=O)OC", True),   # methyl octanoate
        ("CCCCCCCC(=O)OCC(O)CO", True),  # 1-monooctanoylglycerol
        ("CCCCCCCC(=O)O[C@@H](CC([O-])=O)C[N+](C)(C)C", True),  # O-octanoyl-D-carnitine
        ("CCCCCCC(=O)OCC", False),  # ethyl heptanoate (too short)
        ("CCCCCCCCC(=O)OCC", False),  # ethyl nonanoate (too long)
        ("CCCCCCCCOC", False),  # octyl ether (not an ester)
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_octanoate_ester(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()