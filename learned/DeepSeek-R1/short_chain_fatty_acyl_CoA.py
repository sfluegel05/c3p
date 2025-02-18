"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:156875 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA has a CoA moiety linked via thioester to a fatty acid with 2-6 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Improved adenine pattern (purine core with amino group)
    # Matches the adenine component in CoA (N6-amino group)
    adenine_pattern = Chem.MolFromSmarts("[nH0]1c2ncnc2c([nH0])n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety (CoA component)"
    
    # Verify CoA's pantetheine-phosphate backbone components
    # Check for presence of phosphate groups (at least two)
    phosphate_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)  # Phosphorus atoms
    if phosphate_count < 3:  # CoA typically has 3 phosphate groups
        return False, "Insufficient phosphate groups for CoA structure"
    
    # Check for thioester group (C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # Check each thioester group for acyl chain length
    for match in thioester_matches:
        # Match indices: [C, O, S]
        c_carbon = mol.GetAtomWithIdx(match[0])
        
        # Find the R group (non-O/S neighbor of C)
        r_atom = None
        for neighbor in c_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() not in {8, 16}:  # Not O or S
                r_atom = neighbor
                break
        if not r_atom:
            continue  # No R group attached
        
        # Recursively count carbons in R group (excluding carbonyl)
        visited = set()
        def count_carbons(atom):
            if atom.GetIdx() in visited or atom.GetAtomicNum() != 6:
                return 0
            visited.add(atom.GetIdx())
            return 1 + sum(count_carbons(n) for n in atom.GetNeighbors() if n != c_carbon)
        
        r_carbons = count_carbons(r_atom)
        acyl_length = r_carbons + 1  # Include carbonyl carbon
        
        if 2 <= acyl_length <= 6:
            return True, f"Short-chain acyl group ({acyl_length} carbons) attached via thioester to CoA"
    
    return False, "Acyl chain length not in 2-6 carbon range"