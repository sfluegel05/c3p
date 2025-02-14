"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an amino acid where the amine group is acetylated
    (CH3-C(=O)-N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS pattern for N-acetyl amino acid
    acetyl_amino_acid_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3][CX4][CX4]([!H])[CX3](=[OX1])[OX1,OX2-]")
    
    # Check for the presence of the overall substructure
    if not mol.HasSubstructMatch(acetyl_amino_acid_pattern):
      return False, "Molecule does not match the N-acetyl-amino acid core pattern"
      
    # Check for the number of acetyl groups to be one and only one
    acetyl_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3]")
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    
    if len(acetyl_matches) != 1:
      return False, f"Molecule has {len(acetyl_matches)} acetyl groups, must have only 1"

    return True, "Molecule matches the N-acetyl-amino acid pattern."