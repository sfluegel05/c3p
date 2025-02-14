"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    It should match a variety of structures based on the phenylpropanoid class.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broadened search for phenylpropanoid-like structures
    # more generic aromatic system (includes benzene and potential substitutions)
    aromatic_system_pattern = Chem.MolFromSmarts("a1aaaaa1")
    if not mol.HasSubstructMatch(aromatic_system_pattern):
        return False, "No aromatic system found"

    # Try to detect characteristic linkages including various propanoid structures
    phenylpropanoid_broad_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(cc1)-C")  # for C6-C3 backbone clues
    possible_backbones = [
        Chem.MolFromSmarts("c1ccccc1"),  # benzene
        Chem.MolFromSmarts("c1ccc(O)cc1"),  # phenol
        Chem.MolFromSmarts("c1ccccc1C(C)(C)C"),  # cases involving isoprenoids
        Chem.MolFromSmarts("c1ccccc1C=C"),  # phenyl vinyl structures
        phenylpropanoid_broad_pattern
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in possible_backbones):
        return False, "No suitable phenylpropanoid linkage found"

    # Look for functional groups frequently seen in phenylpropanoids
    functional_groups = [
        Chem.MolFromSmarts("[OH]"),  # hydroxyl
        Chem.MolFromSmarts("CO"),  # methoxy group
        Chem.MolFromSmarts("OC=O"),  # ester linkage
        Chem.MolFromSmarts("O=C"),  # carbonyl group
        Chem.MolFromSmarts("C=C"),  # olefinic group
    ]

    has_functional_group = any(mol.HasSubstructMatch(fg) for fg in functional_groups)
    
    if not has_functional_group:
        return False, "Missing characteristic functional groups"

    return True, "Contains aromatic system and phenylpropanoid-associated linkages and functional groups"