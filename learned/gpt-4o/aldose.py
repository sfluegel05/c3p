"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str) -> (bool, str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Aldehyde pattern with at least two adjacent hydroxyl groups
    open_chain_pattern = Chem.MolFromSmarts("C(=O)[CX4H1][CX4H1](O)[CX4H1](O)")

    # Patterns for five- and six-membered cyclic hemiacetals
    furanose_pattern = Chem.MolFromSmarts("O1[C@@H]([CX4H1](O)[CX4H1](O)C1)")
    pyranose_pattern = Chem.MolFromSmarts("O1[C@@H]([CX4H1](O)[C@H](O)[C@H](O)[C@H]1O)")

    # Check open-chain form
    if mol.HasSubstructMatch(open_chain_pattern):
        return True, "Structure is consistent with open-chain form of an aldose"

    # Check five- or six-membered cyclic hemiacetal forms
    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        return True, "Structure is consistent with cyclic hemiacetal form of an aldose"

    return False, "Structure does not fit criteria for an aldose"