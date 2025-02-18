"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine (CHEBI:53438)
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an N-acyl-amino acid where the amino acid is glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains N-acylglycine substructure, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for N-acylglycine core:
    # Acyl group (C=O) attached to nitrogen with one H (excludes peptides),
    # connected to glycine's CH2 and carboxylic acid (protonated or deprotonated)
    pattern = Chem.MolFromSmarts('[CX3](=O)-[NX3H1]-[CH2]-[CX3](=O)-[OX2H1,OX1-]')
    
    if mol.HasSubstructMatch(pattern):
        return True, "Contains N-acylglycine core structure"
    
    # Additional check for esterified carboxylic acid (like in Mortivinacin B)
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)-[NX3H1]-[CH2]-[CX3](=O)-[OX2]')
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Contains N-acylglycine ester variant"
    
    return False, "Does not match N-acylglycine pattern"