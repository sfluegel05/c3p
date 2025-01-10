"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is defined as having an alpha-carbon connected to an amino group
    and a carboxyl group esterified with an alcohol.

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

    # Pattern: Extended alpha amino acid ester structure
    # Including any bond after ester to capture esterified alcohol rather than focusing only on terminal O
    alpha_amino_pattern = "[NX3][C@,C](Cc1~c*~c*1|[R1])"  # Alpha carbon binds an amino group maybe with stereocenter
    ester_linkage_pattern = "[C](=O)O[C,H]"  # Capturing ester linkage with any bound after O
    
    # Combined pattern
    complete_pattern = Chem.MolFromSmarts(alpha_amino_pattern + ester_linkage_pattern)

    if complete_pattern is None:
        return None, None  # Error in pattern creation

    if mol.HasSubstructMatch(complete_pattern):
        return True, "Contains an alpha-amino acid backbone with ester linkage"
    
    return False, "Does not match the alpha-amino acid ester pattern"