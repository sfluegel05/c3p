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
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H,OX1-]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbons and check minimum chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty acid"
        
    # Find different types of branching points
    # Regular carbon branching (sp3)
    branched_carbon_pattern = Chem.MolFromSmarts("[CH1,CH0](-[#6])(-[#6])(-[#6])")
    # Methyl branches
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][CH1]([#6])[#6]")
    # Ethyl or longer branches
    alkyl_branch_pattern = Chem.MolFromSmarts("[CH2][CH1]([#6])[#6]")
    
    branching_points = mol.GetSubstructMatches(branched_carbon_pattern)
    methyl_branches = mol.GetSubstructMatches(methyl_branch_pattern)
    alkyl_branches = mol.GetSubstructMatches(alkyl_branch_pattern)
    
    # Look for cyclopropyl groups (common in mycolic acids)
    cyclopropyl_pattern = Chem.MolFromSmarts("[CH2R1]1[CH1R1][CH1R1]1")
    cyclopropyl_groups = mol.GetSubstructMatches(cyclopropyl_pattern)
    
    if not (branching_points or methyl_branches or alkyl_branches or cyclopropyl_groups):
        return False, "No branching found in carbon chain"
    
    # Check if molecule is primarily aliphatic
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > c_count/2:
        return False, "Molecule is primarily aromatic"
    
    # Look for a substantial carbon chain
    # Count rotatable bonds as proxy for chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Carbon chain too rigid/short for fatty acid"
    
    # Additional checks for reasonable fatty acid properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 60:  # Minimum weight for a small branched fatty acid
        return False, "Molecular weight too low for branched-chain fatty acid"
    
    # Count oxygen atoms (should typically be 2 or more)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for fatty acid"
        
    # Build reason string based on structural features
    features = []
    if methyl_branches:
        features.append(f"{len(methyl_branches)} methyl branch(es)")
    if alkyl_branches:
        features.append(f"{len(alkyl_branches)} alkyl branch(es)")
    if cyclopropyl_groups:
        features.append(f"{len(cyclopropyl_groups)} cyclopropyl group(s)")
    
    reason = "Contains carboxylic acid group and " + ", ".join(features)
    if o_count > 2:
        reason += f" and {o_count-2} additional oxygen(s)"
    
    return True, reason