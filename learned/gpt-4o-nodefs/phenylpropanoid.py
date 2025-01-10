"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is likely a phenylpropanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for phenyl group (benzene ring)
    phenyl_group = Chem.MolFromSmarts('c1ccccc1')
    if not mol.HasSubstructMatch(phenyl_group):
        return False, "No phenyl group (benzene ring) found"

    # Check for various propanoid derivatives (common SMILES patterns observed in phenylpropanoids)
    phenylpropanoid_like_patterns = [
        Chem.MolFromSmarts('c1ccc(cc1)C(=O)O'),  # Cinnamic acid derivatives
        Chem.MolFromSmarts('c1ccc(cc1)C(O)C=O'), # Coniferyl aldehyde
        Chem.MolFromSmarts('c1ccc(cc1)C(O)C[OH]'), # Alcohols derived from coniferyl alcohol
        Chem.MolFromSmarts('c1ccc(cc1)COC=O'),  # Coumarins
        Chem.MolFromSmarts('c1ccc(cc1)C=C[CH2]'),  # Ferulic acid
    ]

    found_phenylpropanoid_like_pattern = any(mol.HasSubstructMatch(pattern) for pattern in phenylpropanoid_like_patterns)
    if not found_phenylpropanoid_like_pattern:
        return False, "No recognizable phenylpropanoid-like pattern found"

    # Check for common functional groups such as hydroxyls, methoxy, carbonyl, etc.
    functional_groups = [
        Chem.MolFromSmarts('[OH]'),   # Hydroxyl
        Chem.MolFromSmarts('COC'),    # Methoxy
        Chem.MolFromSmarts('[CX3]=[O]'),  # Carbonyl in ketones/aldehydes
        Chem.MolFromSmarts('C(=O)O'),  # Ester linkages
        Chem.MolFromSmarts('CO'),     # Ether linkages
    ]

    if not any(mol.HasSubstructMatch(fg) for fg in functional_groups):
        return False, "No characteristic functional groups found"

    return True, "Contains phenyl group and derivative patterns frequently found in phenylpropanoids"