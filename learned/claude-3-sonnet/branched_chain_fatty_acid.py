"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:35741 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbons and check minimum chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty acid"
        
    # Find branching points (carbons with >2 carbon neighbors)
    branched_carbon_pattern = Chem.MolFromSmarts("[CH1,CH0](-[#6])(-[#6])(-[#6])")
    branching_points = mol.GetSubstructMatches(branched_carbon_pattern)
    
    # Also look for methyl branches specifically
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][CH1]([#6])[#6]")
    methyl_branches = mol.GetSubstructMatches(methyl_branch_pattern)
    
    if not (branching_points or methyl_branches):
        return False, "No branching found in carbon chain"
    
    # Check if molecule is primarily aliphatic (exclude primarily aromatic compounds)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > c_count/2:
        return False, "Molecule is primarily aromatic"
    
    # Additional checks for reasonable fatty acid properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 60:  # Minimum weight for a small branched fatty acid
        return False, "Molecular weight too low for branched-chain fatty acid"
    
    # Count oxygen atoms (should be 2 for standard fatty acid)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for fatty acid"
        
    # Look for a continuous carbon chain
    longest_chain = rdMolDescriptors.CalcLongestChain(mol)
    if longest_chain < 4:
        return False, "Main carbon chain too short"
    
    # Success case - specify the type of branching found
    reason = "Contains carboxylic acid group and "
    if branching_points and methyl_branches:
        reason += f"multiple branch points including {len(methyl_branches)} methyl branches"
    elif branching_points:
        reason += f"{len(branching_points)} branch points"
    else:
        reason += f"{len(methyl_branches)} methyl branches"
    
    return True, reason