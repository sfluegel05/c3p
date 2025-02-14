"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an alpha-amino acid core where the carboxyl group
    is replaced by an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for alpha-amino acid ester core.
    # [C@H](N)([CX3](=O)[OX2][CX4]) matches a chiral carbon (alpha carbon) attached to:
    #    - an amine (N)
    #    - a carbonyl (C=O), that connects to
    #    - an ester (OX2[CX4])
    # This is more specific than before because it is looking for connectivity.
    alpha_amino_ester_pattern = Chem.MolFromSmarts("[C@H]([NX3H0-2])([CX3](=[OX1])[OX2][#6])")


    #check that the pattern is matched
    alpha_amino_matches = mol.GetSubstructMatches(alpha_amino_ester_pattern)


    if not alpha_amino_matches:
        return False, "Molecule does not contain an alpha-amino acid ester core"

    return True, "Molecule contains an alpha-amino acid with an ester group."