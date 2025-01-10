"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should typically have just one
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches > 2:  # allowing up to 2 for some edge cases
        return False, f"Too many carboxylic acid groups ({carboxylic_matches})"
    
    # Look for methyl branches (excluding those part of carboxylic group)
    # Pattern matches methyl group attached to carbon that's not part of COOH
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][CX4]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_branch_pattern))
    if methyl_matches == 0:
        return False, "No methyl branches found"
    
    # Check for non-methyl branches (e.g. ethyl, propyl)
    # This pattern matches longer branches (C-C-C or more)
    long_branch_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2,CH3]")
    branch_matches = mol.GetSubstructMatches(long_branch_pattern)
    
    # Count carbons and check if it's a reasonable size for a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:  # minimum size for a branched fatty acid
        return False, "Carbon chain too short"
    
    # Check for predominantly aliphatic character
    aromatic_pattern = Chem.MolFromSmarts("a")
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings"
    
    # Count heteroatoms (excluding O from COOH)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Allow some nitrogen-containing compounds but limit other heteroatoms
    if s_count > 0 or p_count > 0:
        return False, "Contains unexpected heteroatoms (S, P)"
    if n_count > 2:  # allowing up to 2 nitrogens for some derivatives
        return False, "Too many nitrogen atoms"
    
    # Calculate fraction of carbons that are sp3 hybridized
    sp3_pattern = Chem.MolFromSmarts("[CX4]")
    sp3_count = len(mol.GetSubstructMatches(sp3_pattern))
    sp3_fraction = sp3_count / c_count
    
    # Most carbons should be sp3 (allowing some unsaturation)
    if sp3_fraction < 0.5:
        return False, "Too many unsaturated carbons"
    
    return True, f"Contains carboxylic acid group with {methyl_matches} methyl branch(es)"