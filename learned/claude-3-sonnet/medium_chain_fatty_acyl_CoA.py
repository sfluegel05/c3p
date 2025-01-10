"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA moiety patterns
    # Adenine base
    adenine_pattern = Chem.MolFromSmarts("n1c(N)nc2c(ncnc12)")
    # Phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    # Pantetheine part
    pantetheine_pattern = Chem.MolFromSmarts("CC(C)(COP)C(O)C(=O)NCCC(=O)NCCS")
    
    if not all(mol.HasSubstructMatch(p) for p in [adenine_pattern, phosphate_pattern, pantetheine_pattern]):
        return False, "Missing essential CoA structural elements"
    
    # Check for thioester linkage (R-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Count carbons in the fatty acid chain
    # First, get the thioester carbon
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Could not analyze fatty acid chain"
    
    # Get the carbon atom index of the thioester carbonyl
    carbonyl_idx = thioester_matches[0][0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Count carbons in the chain (excluding the carbonyl)
    chain_carbons = set()
    def count_chain_carbons(atom, visited):
        if atom.GetAtomicNum() != 6:  # not carbon
            return
        if atom.GetIdx() in visited:
            return
        visited.add(atom.GetIdx())
        chain_carbons.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() != carbonyl_idx:  # don't go back through carbonyl
                count_chain_carbons(neighbor, visited)
    
    # Start counting from carbons attached to carbonyl
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # carbon
            count_chain_carbons(neighbor, set())
    
    chain_length = len(chain_carbons)
    
    if chain_length < 6:
        return False, f"Fatty acid chain too short ({chain_length} carbons)"
    if chain_length > 12:
        return False, f"Fatty acid chain too long ({chain_length} carbons)"
        
    # Additional checks for reasonable molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 750 or mol_weight > 1100:
        return False, f"Molecular weight {mol_weight:.1f} outside typical range"
    
    return True, f"Medium-chain fatty acyl-CoA with {chain_length} carbons in fatty acid chain"