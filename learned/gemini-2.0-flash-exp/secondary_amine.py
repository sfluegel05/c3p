"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to two carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary amine
    # [NX3;!H0]([C])([C]) specifies a nitrogen with 3 connections (including implicit hydrogens), at least one hydrogen (not zero) with two carbons.
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;!H0]([C])([C])")
    

    # Check for at least one match
    if not mol.HasSubstructMatch(secondary_amine_pattern):
        return False, "No secondary amine found (N with two carbons)"


    return True, "Molecule contains a secondary amine group"