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

    # Define SMARTS pattern for N-acetyl-amino acid
    # The pattern searches for: Acetyl group (CC(=O)N), connected to an alpha-carbon (C), part of amino acid (C(N)C(=O)O)
    acetyl_amino_acid_pattern = Chem.MolFromSmarts("C(C(=O)O)[NH]C(=O)C")
    
    # Check for the presence of the N-acetyl-amino acid pattern
    if mol.HasSubstructMatch(acetyl_amino_acid_pattern):
        return True, "Contains N-acetyl group attached to nitrogen in an amino acid structure"
    else:
        return False, "Does not contain the characteristic N-acetyl group in amino acid configuration"