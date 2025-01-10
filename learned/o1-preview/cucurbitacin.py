"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.
    They have a characteristic cucurbitane skeleton with specific stereochemistry.

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

    # Define the cucurbitane skeleton SMARTS pattern with stereochemistry
    cucurbitane_smarts = """
    [#6@H]1CC[C@@]2(C)[C@@H](CC[C@]3(C)[C@]2(CC[C@@H]1C)C)[C@H]1CC[C@@]4(C)[C@@H](CC[C@]4(C)[C@]3(CC1)C)C
    """
    cucurbitane_pattern = Chem.MolFromSmarts(cucurbitane_smarts)
    if cucurbitane_pattern is None:
        return False, "Invalid cucurbitane SMARTS pattern"

    # Check for cucurbitane skeleton
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "Does not contain the cucurbitane skeleton characteristic of cucurbitacins"

    # Check for characteristic functional groups
    # Cucurbitacins often have multiple hydroxyl and ketone groups, and an alpha,beta-unsaturated ketone
    num_hydroxyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C;H1,H2]-[OH]')))
    num_ketone = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[C;!$(C=O)]')))
    num_enone = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C=CC=O')))

    if num_hydroxyl + num_ketone + num_enone < 3:
        return False, f"Contains {num_hydroxyl} hydroxyls, {num_ketone} ketones, and {num_enone} enones, less than 3 total characteristic functional groups"

    # Check for triterpenoid skeleton (30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Carbon count {c_count} is less than 30, not a triterpenoid"

    # Additional check for common functional groups in cucurbitacins
    # Cucurbitacins often have an acetoxy group at C25
    acetoxy_pattern = Chem.MolFromSmarts('C(C)(C)OC(=O)C')
    if len(mol.GetSubstructMatches(acetoxy_pattern)) > 0:
        return True, "Contains cucurbitane skeleton and characteristic functional groups of cucurbitacins, including acetoxy group"

    return True, "Contains cucurbitane skeleton and characteristic functional groups of cucurbitacins"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cucurbitacin',
        'definition': 'Any one of a class of tetracyclic triterpenoids, formally derived from the triterpene hydrocarbon cucurbitane, developed by some plants (especially those of the family Cucurbitaceae) as a defence mechanism against herbivores.',
        'parents': []
    }
}