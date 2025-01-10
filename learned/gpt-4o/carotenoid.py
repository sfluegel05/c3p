"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are C40 tetraterpenoids with a long conjugated double-bond system,
    potentially containing cyclic structures, oxygen functionalities, and excluding retinoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 40:
        return False, f"Number of carbon atoms is {c_count}, should be 40 for carotenoids"

    # Check for presence of conjugated double bonds
    conjugated_double_bond_pattern = Chem.MolFromSmarts("C=C")
    conjugated_double_bond_matches = mol.GetSubstructMatches(conjugated_double_bond_pattern)
    if len(conjugated_double_bond_matches) < 8:
        return False, f"Fewer conjugated double bonds found ({len(conjugated_double_bond_matches)}), suggestive of non-carotenoid structure"

    # Check for cyclic structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No cyclic structures found, which are common in carotenoids"

    # Exclude retinoid-like structures (simplified pattern match as example)
    retinoid_pattern = Chem.MolFromSmarts("CC1=C(CCCC1)C(C)=C")
    if mol.HasSubstructMatch(retinoid_pattern):
        return False, "Structure matches common retinoid pattern, excluded from carotenoids"

    # Check for presence of functional groups (for xanthophylls)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 0:
        return True, "Classified as a xanthophyll carotenoid due to the presence of oxygen functionality"
        
    return True, "Classified as a carotene (oxygen-less carotenoid) due to C40 tetraterpenoid structure with conjugated bonds"