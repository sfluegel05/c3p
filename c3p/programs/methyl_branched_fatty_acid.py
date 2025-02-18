"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: CHEBI:35839 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid with methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for atoms other than C, H, and O
    allowed_atoms = [6, 1, 8]  # C, H, O
    if any(atom.GetAtomicNum() not in allowed_atoms for atom in mol.GetAtoms()):
        return False, "Contains atoms other than C, H, and O"
    
    # Check for methyl branches
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_branch_matches = mol.GetSubstructMatches(methyl_branch_pattern)
    
    # Check if methyl branches are attached to non-terminal carbons
    methyl_branches_on_chain = False
    for match in methyl_branch_matches:
        methyl_carbon = mol.GetAtomWithIdx(match[0])
        neighbors = [mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() for neighbor_idx in methyl_carbon.GetNeighbors()]
        if 6 in neighbors:  # Check for carbon neighbor
            methyl_branches_on_chain = True
            break
    
    if not methyl_branches_on_chain:
        return False, "No methyl branches attached to non-terminal carbons"
    
    # Check for other types of branches
    other_branch_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]")
    other_branch_matches = mol.GetSubstructMatches(other_branch_pattern)
    if other_branch_matches:
        return False, "Contains branches other than methyl groups"
    
    return True, "Molecule is a methyl-branched fatty acid"