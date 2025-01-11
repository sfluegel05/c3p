"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is defined by having an acetyl group as the acyl group attached to the nitrogen.

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
    # Looking for: Acetyl group (C(=O)C) attached to a nitrogen (N)
    acetyl_amino_pattern = Chem.MolFromSmarts("CC(=O)N")
    
    # Look for additional carboxylate group (C(O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Consider carboxylate ion form
    carboxylic_acid_neutral_pattern = Chem.MolFromSmarts("C(=O)O")  # Consider neutral form
    
    if (mol.HasSubstructMatch(acetyl_amino_pattern) and 
        (mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(carboxylic_acid_neutral_pattern))):
        return True, "Contains acetyl group attached to nitrogen with carboxylic acid group, classifies as N-acetyl-amino acid"
    else:
        return False, "Does not contain required acetyl-to-nitrogen bond with adjacent carboxylic acid group"