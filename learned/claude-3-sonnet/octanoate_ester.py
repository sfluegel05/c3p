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

    # Look for ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Pattern for octanoate group:
    # - 8 carbons in chain (including carbonyl carbon)
    # - Ends in ester group
    # - More flexible pattern that allows for substitutions and stereochemistry
    octanoate_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX3](=[OX1])[OX2]")
    
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No octanoate chain found"
    
    matches = mol.GetSubstructMatches(octanoate_pattern)
    
    for match in matches:
        # Get the atoms in the chain
        chain_atoms = list(match[:-2])  # Exclude ester oxygen and carbonyl oxygen
        
        # Check that the carbons in the chain form a continuous path
        is_continuous = True
        for i in range(len(chain_atoms)-1):
            bond = mol.GetBondBetweenAtoms(chain_atoms[i], chain_atoms[i+1])
            if bond is None:
                is_continuous = False
                break
                
        if not is_continuous:
            continue
            
        # Check that the chain itself isn't part of a ring
        chain_in_ring = False
        for i in range(len(chain_atoms)-1):
            bond = mol.GetBondBetweenAtoms(chain_atoms[i], chain_atoms[i+1])
            if bond.IsInRing():
                chain_in_ring = True
                break
                
        if chain_in_ring:
            continue
            
        # Verify the chain length (should be 8 carbons including carbonyl)
        chain_length = len(chain_atoms) + 1  # Add 1 for carbonyl carbon
        if chain_length != 8:
            continue
            
        # Check that all atoms in chain are carbons
        all_carbons = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in chain_atoms)
        if not all_carbons:
            continue
            
        # If we get here, we've found a valid octanoate ester group
        return True, "Contains octanoate ester group (8-carbon chain with terminal ester)"
    
    return False, "No valid octanoate ester group found"

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