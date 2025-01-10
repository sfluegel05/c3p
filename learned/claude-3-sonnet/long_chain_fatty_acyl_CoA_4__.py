"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA core structure - adenine base
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety of CoA"
    
    # Check for thioester linkage (R-C(=O)-S-R)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatches(thioester_pattern):
        return False, "Missing thioester linkage"
    
    # Check for phosphate groups - multiple patterns to catch different representations
    phosphate_patterns = [
        Chem.MolFromSmarts("P(=O)([O-])[O-]"),
        Chem.MolFromSmarts("P(=O)(O)O"),
        Chem.MolFromSmarts("P(=O)(O)[O-]")
    ]
    phosphate_count = 0
    for pattern in phosphate_patterns:
        phosphate_count += len(mol.GetSubstructMatches(pattern))
    if phosphate_count < 3:
        return False, f"Found only {phosphate_count} phosphate groups, need at least 3"
    
    # Check for pantetheine portion
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine portion of CoA"
    
    # Look for long carbon chain - multiple patterns to catch different types
    chain_patterns = [
        # Saturated chain
        Chem.MolFromSmarts("CCCCCCCCCCC"),
        # Chain with possible double bonds
        Chem.MolFromSmarts("C~C~C~C~C~C~C~C~C~C~C"),
        # Chain with possible branches
        Chem.MolFromSmarts("[CH2,CH3]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]")
    ]
    
    has_long_chain = False
    for pattern in chain_patterns:
        if mol.HasSubstructMatch(pattern):
            has_long_chain = True
            break
            
    if not has_long_chain:
        return False, "No long carbon chain found (need at least 11 carbons)"
    
    # Identify features
    features = []
    
    # Check for unsaturation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("unsaturated")
        
    # Check for hydroxyl groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2,CH]-[OH]")):
        features.append("hydroxylated")
        
    # Check for keto groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)C")):
        features.append("keto")
        
    # Check for epoxy groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        features.append("epoxy")
    
    feature_str = " and ".join(features) if features else "saturated"
    
    # Calculate approximate chain length
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        # Get the carbon attached to the thioester
        thioester_carbon = thioester_matches[0][0]
        # Count connected carbons in the chain
        visited = set()
        def count_chain_carbons(atom_idx, depth=0):
            if depth > 30 or atom_idx in visited:  # Prevent infinite recursion
                return 0
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:  # Not carbon
                return 0
            count = 1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    count += count_chain_carbons(neighbor.GetIdx(), depth + 1)
            return count
            
        chain_length = count_chain_carbons(thioester_carbon)
    else:
        chain_length = "unknown"

    return True, f"Long-chain ({chain_length} carbons) {feature_str} fatty acyl-CoA(4-) with thioester linkage"