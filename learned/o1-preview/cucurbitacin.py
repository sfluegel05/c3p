"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""

from rdkit import Chem

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.
    They have a characteristic tetracyclic skeleton with specific functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove stereochemistry from the molecule to make matching less specific
    Chem.RemoveStereochemistry(mol)

    # Use the cucurbitane backbone as substructure pattern
    # Cucurbitane SMILES without stereochemistry
    cucurbitane_smiles = 'CC1CCC2(C)C1CCC3(C)C(C)CCC4C(C)(C)CCC23C4'
    cucurbitane_mol = Chem.MolFromSmiles(cucurbitane_smiles)
    if cucurbitane_mol is None:
        return False, "Invalid cucurbitane SMILES"

    # Remove stereochemistry from the cucurbitane molecule
    Chem.RemoveStereochemistry(cucurbitane_mol)

    # Check for cucurbitane skeleton
    if not mol.HasSubstructMatch(cucurbitane_mol):
        return False, "Does not contain the cucurbitane skeleton characteristic of cucurbitacins"

    # Check for characteristic functional groups
    # Cucurbitacins often have multiple hydroxyl groups and ketones
    num_hydroxyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]')))
    num_ketone = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[!$([O-])]')))

    total_functional_groups = num_hydroxyl + num_ketone
    if total_functional_groups < 3:
        return False, f"Contains {total_functional_groups} hydroxyls and ketones, less than 3 total characteristic functional groups"

    # Check for triterpenoid skeleton (at least 30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Carbon count {c_count} is less than 30, not a triterpenoid"

    return True, "Contains cucurbitane skeleton and characteristic functional groups of cucurbitacins"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cucurbitacin',
        'definition': 'Any one of a class of tetracyclic triterpenoids, formally derived from the triterpene hydrocarbon cucurbitane, developed by some plants (especially those of the family Cucurbitaceae) as a defence mechanism against herbivores.',
        'parents': []
    }
}