"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid contains a steroid core structure with a hydroxyl group at the 11beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved steroid core pattern recognition: wide steric variety including stereocenters
    steroid_core_patterns = [
        Chem.MolFromSmarts("C[C@@H]1CC[C@@H]2[C@H]3CC[C@H]4CCC(=O)C=C[C@@]4(C)[C@@H]3CCC2C1"),  # Common stereo arrangement
        Chem.MolFromSmarts("C1[C@H]2[C@H](C[C@H]3[C@H]2CC[C@@H]4C3(C)CC=C4)CC[C@@H]1C=O"),  # Alternate stereo
        Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C4)C3CCC2C1")  # General steroid scaffold
    ]

    found_steroid = any(mol.HasSubstructMatch(p) for p in steroid_core_patterns)
    if not found_steroid:
        return False, "Steroid core structure not found"

    # Specific 11beta-hydroxy group: looking for hydroxyl on specific positions
    hydroxy_patterns = [
        Chem.MolFromSmarts("CC1(C)C[C@H](O)[C@@H]2CCC3C(CCC4=[O])C3C[C@@H]2C(=C1)O"),  # Comprehensive stereochem hydroxyl specification
        Chem.MolFromSmarts("C[C@@H]1C[C@@H](O)C2Cc3cc(O)ccc3C2CC1"),  # Less specific steric focus
        Chem.MolFromSmarts("[C@@H]([C@H](O)[C@@H]2CCC3C4C(CC[C@]34C)[C@H]2[C@@H](C)C1)O")  # Complete stereo-specific signature
    ]

    found_hydroxy = any(mol.HasSubstructMatch(p) for p in hydroxy_patterns)
    if not found_hydroxy:
        return False, "Missing 11beta-hydroxy group"

    return True, "Contains steroid core structure with 11beta-hydroxy group"