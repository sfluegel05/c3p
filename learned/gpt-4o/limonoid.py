"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid with a structure commonly
    derived from a 4,4,8-trimethyl-17-furanylsteroid skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a reasonable range of carbon atoms typical in limonoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 40:
        return False, f"Number of carbon atoms {c_count} outside expected range for limonoid"

    # Ensure a high range of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Insufficient oxygen content for a highly oxygenated limonoid: {o_count}"

    # Identify oxygen functional groups: ketones, esters, hydroxyl, and potential lactone
    important_oxy_patterns = [
        Chem.MolFromSmarts("C=O"),
        Chem.MolFromSmarts("C(=O)O"),
        Chem.MolFromSmarts("[OX2H]"),  # Alcohol
        Chem.MolFromSmarts("O1CCO1"),  # Broader to catch different heterocycles
        Chem.MolFromSmarts("C1(O)CC1"),  # Another form of closed-ring oxygen groups
    ]
    for oxy_pattern in important_oxy_patterns:
        if mol.HasSubstructMatch(oxy_pattern):
            break
    else:
        return False, "Missing characteristic oxygenated functionalities"

    # Check if there is any form of methyl groups, allowing for variations
    methyl_patterns = [
        Chem.MolFromSmarts("C(C)(C)C"),  # Standard trimethyl
        Chem.MolFromSmarts("C(C)(C)(C)"),  # Various methylated forms
        Chem.MolFromSmarts("C(C)([OX2])")  # Methyl ether or alcohol derivatives
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in methyl_patterns):
        return False, "Lack of characteristic methylation"

    return True, "SMILES string corresponds to a limonoid structure"