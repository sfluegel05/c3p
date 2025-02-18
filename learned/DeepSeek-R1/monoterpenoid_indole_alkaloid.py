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
    
    # Improved indole detection: fused benzene + pyrrole (allowing substitutions)
    indole_pattern = Chem.MolFromSmarts("[#7]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~2~1")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety detected"
    
    # Check for monoterpenoid features (isoprenoid patterns)
    # Look for at least two isoprene units (C10 skeleton) - approximate via branching
    # Pattern for head-to-tail isoprene units: CH2-C(CH2)-CH2-CH2 (approximate)
    isoprene_pattern = Chem.MolFromSmarts("[CH2][CH]([CH2])[CH2][CH2]")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No diisoprenoid (terpenoid) features detected"
    
    # Check for typical linkage between indole and terpenoid (e.g., ester/amine)
    ester_amine_linkage = Chem.MolFromSmarts("[#7,#8]-[#6]-[#6](-[#6](-[#6])-[#6])")  # Approx
    if not mol.HasSubstructMatch(ester_amine_linkage):
        return False, "Missing indole-terpenoid linkage"
    
    # Adjusted carbon count (monoterpene C10 + indole C9 ~19, allow some variation)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:  # Some examples have lower (e.g., yohimban has 19)
        return False, f"Insufficient carbons ({c_count}) for monoterpenoid indole alkaloid"
    
    # Check for nitrogen presence (alkaloid characteristic)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen atoms (not an alkaloid)"
    
    return True, "Contains indole moiety with diisoprenoid features and alkaloid characteristics"