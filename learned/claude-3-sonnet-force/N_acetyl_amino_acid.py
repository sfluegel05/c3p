"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:38756 N-acetyl-amino acid
An N-acyl-amino acid that has acetyl as the acyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.

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
    
    # Look for acetyl group attached to nitrogen
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group attached to nitrogen"
    
    # Look for amino group (attached to the nitrogen with the acetyl group)
    amino_pattern = Chem.MolFromSmarts("[N;X3][C;X4]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    
    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check if the nitrogen attached to the acetyl group is also part of the amino group
    acetyl_N = mol.GetSubstructMatches(acetyl_pattern)[0][1]
    amino_N = mol.GetSubstructMatches(amino_pattern)[0][1]
    if acetyl_N != amino_N:
        return False, "Acetyl group not attached to amino nitrogen"
    
    return True, "Contains acetyl group attached to nitrogen of an amino acid"