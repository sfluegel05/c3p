"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

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

    # Look for carboxylic acid group -COOH (smarts pattern: C(=O)[O;H])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
        
    # Look for methyl branches (smarts pattern: [CH3X4] - methyl group)
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3X4]")
    methyl_branch_matches = mol.GetSubstructMatches(methyl_branch_pattern)
    if len(methyl_branch_matches) < 1:
        return False, "No methyl branches found"
    
    # Count number of carbon atoms to confirm it's a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6: # Fatty acids typically have at least 6 carbons
        return False, "Not enough carbon atoms (need at least 6 for a fatty acid)"
    
    return True, "Contains carboxylic acid group and has methyl branches"