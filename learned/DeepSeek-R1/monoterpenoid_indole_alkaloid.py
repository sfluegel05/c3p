"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: CHEBI:65312 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Must contain an indole moiety and structural features from a monoterpenoid (diisoprenoid).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Corrected indole detection: fused benzene and pyrrole (nH in 5-membered ring)
    indole_pattern = Chem.MolFromSmarts("[nH]1ccc2ccccc12")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety detected"
    
    # Check for monoterpenoid features (approximate C10 skeleton with branching)
    # Look for at least two methyl groups attached to non-terminal carbons as proxy for isoprenoid
    methyl_branch = Chem.MolFromSmarts("[CH3;!$(C[OH,O,NH1,NH2])]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_branch))
    if methyl_matches < 2:
        return False, "Insufficient methyl groups for monoterpenoid"
    
    # Check for typical linkage between indole and terpenoid (e.g., C-C bridge)
    # Pattern matches a chain of at least 3 carbons connecting two rings
    linker_pattern = Chem.MolFromSmarts("[*]!@[CX4]!@[CX4]!@[CX4]!@[*]")
    if not mol.HasSubstructMatch(linker_pattern):
        return False, "Missing indole-terpenoid linkage"
    
    # Adjusted carbon count (monoterpene C10 + indole C9 ~19, allow higher)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 19:
        return False, f"Insufficient carbons ({c_count}) for monoterpenoid indole alkaloid"
    
    # Check for nitrogen presence (alkaloid characteristic)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen atoms (not an alkaloid)"
    
    return True, "Contains indole moiety with monoterpenoid features and alkaloid characteristics"