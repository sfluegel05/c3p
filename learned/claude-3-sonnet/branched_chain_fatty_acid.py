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

    # Check for peptide bonds - exclude peptides
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide bonds"
        
    # Count heteroatoms (excluding O from carboxylic acid)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if n_count > 1 or s_count > 0:
        return False, "Too many heteroatoms for fatty acid"
    
    # Find different types of branching points
    # Regular carbon branching (sp3)
    branched_carbon_pattern = Chem.MolFromSmarts("[CH1,CH0](-[#6])(-[#6])(-[#6])")
    # Methyl branches (including iso- and anteiso-)
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][CH1]([#6])[#6]")
    # Terminal dimethyl branching
    terminal_dimethyl = Chem.MolFromSmarts("[CH3][C]([CH3])([#6])[#6]")
    # Unsaturated branches
    unsaturated_branch = Chem.MolFromSmarts("[CH2,CH3][CH1]=[C]")
    
    branching_points = mol.GetSubstructMatches(branched_carbon_pattern)
    methyl_branches = mol.GetSubstructMatches(methyl_branch_pattern)
    terminal_branches = mol.GetSubstructMatches(terminal_dimethyl)
    unsaturated = mol.GetSubstructMatches(unsaturated_branch)
    
    # Look for cyclopropyl groups (common in mycolic acids)
    cyclopropyl_pattern = Chem.MolFromSmarts("[CH2R1]1[CH1R1][CH1R1]1")
    cyclopropyl_groups = mol.GetSubstructMatches(cyclopropyl_pattern)
    
    # Check if any type of branching is present
    has_branching = (branching_points or methyl_branches or terminal_branches or 
                    cyclopropyl_groups or unsaturated)
    if not has_branching:
        return False, "No branching found in carbon chain"
    
    # Check if molecule is primarily aliphatic
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 2:  # Allow small aromatic substituents
        return False, "Molecule is too aromatic for fatty acid"
    
    # Count rings - fatty acids should have few rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2 and not cyclopropyl_groups:
        return False, "Too many rings for fatty acid"
    
    # Additional checks for reasonable fatty acid properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 50:  # Lower minimum weight to catch smaller acids
        return False, "Molecular weight too low for branched-chain fatty acid"
    
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for fatty acid"
        
    # Build reason string based on structural features
    features = []
    if methyl_branches or terminal_branches:
        n_methyl = len(methyl_branches) + 2*len(terminal_branches)
        features.append(f"{n_methyl} methyl branch(es)")
    if branching_points:
        features.append(f"{len(branching_points)} branching point(s)")
    if cyclopropyl_groups:
        features.append(f"{len(cyclopropyl_groups)} cyclopropyl group(s)")
    if unsaturated:
        features.append(f"{len(unsaturated)} unsaturated branch(es)")
    
    reason = "Contains carboxylic acid group"
    if features:
        reason += " with " + ", ".join(features)
    if o_count > 2:
        reason += f" and {o_count-2} additional oxygen(s)"
    
    return True, reason