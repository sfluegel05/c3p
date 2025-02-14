"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is a glycine where the nitrogen atom is acylated with another group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acylglycine, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Corrected N-acylglycine pattern
    acyl_glycine_pattern = Chem.MolFromSmarts("[N]-[C](=[O])-[C](-[H])(-[H])-[C](=O)[O]")

    if not mol.HasSubstructMatch(acyl_glycine_pattern):
        return False, "Molecule does not contain the N-acylglycine core structure"

    return True, "Molecule contains a glycine with the nitrogen acylated (N-acylglycine)"