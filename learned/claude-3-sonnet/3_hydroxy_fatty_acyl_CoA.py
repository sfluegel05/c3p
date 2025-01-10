"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA core structure
    coenzyme_a = Chem.MolFromSmarts("[#6]-3-[#6](-[#6](-[#6]-[#8]-[P;X4](-[#8])(-[#8])-[#8]-[P;X4](-[#8])(-[#8])-[#8]-[#6]-[#6]-4-[#8]-[#6](-[#6](-[#8])-[#6]-4-[#8]-[P;X4](-[#8])(-[#8])=[#8])-n:2:c:1:c(:n:c(:n:c:1:[#7]):n:c:2)-[#7])(-[#6])-[#6])-[#8])-[#7]-[#6]-[#6]-[#6](-[#7]-[#6]-[#6]-[S]-)=O")
    if not mol.HasSubstructMatch(coenzyme_a):
        return False, "No complete CoA structure found"

    # Look for thioester linkage
    thioester = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.GetSubstructMatches(thioester):
        return False, "No thioester linkage found"

    # Check for 3-hydroxy pattern near thioester with various substitution patterns
    hydroxy_patterns = [
        # Linear chains
        "[CH2,CH3]-[CH1,CH0]([OH1,OH0-])-[CH2X4]-C(=O)S", # Generic
        "[CH2,CH3]-[C@H,C@@H,C]([OH1,OH0-])-[CH2X4]-C(=O)S", # With stereochemistry
        # Branched chains
        "[CH2,CH3,CH]-[C]([OH1,OH0-])([#6])-[CH2X4]-C(=O)S", # Quaternary carbon
        "[CH2,CH3]-[CH1,CH0]([OH1,OH0-])-[CH1,CH0]([#6])-C(=O)S", # Branched at beta position
        # Ring systems
        "[#6]-1-[#6]-[#6]-[#6]([OH1,OH0-])-[#6]-[#6]-1-C(=O)S" # Cyclohexane derivatives
    ]
    
    has_hydroxy = False
    hydroxy_match = None
    for pattern in hydroxy_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if matches:
            has_hydroxy = True
            hydroxy_match = matches[0]
            break
            
    if not has_hydroxy:
        return False, "No 3-hydroxy group found in correct position relative to thioester"

    # Verify the position of the hydroxy group
    if hydroxy_match:
        # Get the carbon bearing the OH group
        hydroxy_carbon = hydroxy_match[1]  # Index may need adjustment based on pattern
        
        # Find path to thioester
        thioester_matches = mol.GetSubstructMatches(thioester)
        if thioester_matches:
            thioester_carbon = thioester_matches[0][0]
            # Check if hydroxy is 3 carbons away from thioester carbon
            path = Chem.GetShortestPath(mol, hydroxy_carbon, thioester_carbon)
            if len(path) != 4:  # Should be 3 bonds apart (4 atoms in path)
                return False, "Hydroxy group not in 3-position relative to thioester"

    # Count carbons in fatty acid portion
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2,CH3]-[CH1,CH0]([OH1,OH0-])-[CH2X4]-C(=O)S")
    if mol.HasSubstructMatch(fatty_acid_pattern):
        matches = mol.GetSubstructMatches(fatty_acid_pattern)
        start_atom = matches[0][0]  # First carbon of pattern
        
        # Count continuous carbon chain from start
        visited = set()
        carbon_count = 0
        
        def count_carbons(atom_idx):
            if atom_idx in visited:
                return
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                nonlocal carbon_count
                carbon_count += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                        count_carbons(neighbor.GetIdx())
        
        count_carbons(start_atom)
        
        if carbon_count < 4:
            return False, f"Carbon chain too short for fatty acid (found {carbon_count} carbons)"

    return True, "Contains CoA moiety with 3-hydroxy fatty acid chain attached via thioester linkage"