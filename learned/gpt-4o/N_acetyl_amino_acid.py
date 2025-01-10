"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is defined by having an acetyl group as the acyl group attached to the nitrogen of an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for acetyl group directly attached to an amino nitrogen
    acetyl_to_amino_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)[O,O-])C(=O)C")
    
    # Check for the presence of the pattern
    if mol.HasSubstructMatch(acetyl_to_amino_pattern):
        return True, "Contains acetyl group directly attached to nitrogen with a typical amino acid structure"
    else:
        return False, "Does not contain the characteristic N-acetyl group in amino acid configuration"