"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group linked to the nitrogen of an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for N-acetyl group pattern (N-C(=O)-C)
    acetyl_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)C")
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No N-acetyl group found"
    
    # Look for amino acid pattern − must have a C(=O)O and N attached to chiral alpha-carbon center
    amino_acid_pattern = Chem.MolFromSmarts("C(=O)(O)[C@H]N")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid structure detected"

    return True, "Contains N-acetyl group attached to an amino acid structure"