"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI: CHEBI:76295 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an N-acyl-amino acid with acetyl as the acyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"
    
    # Check for acetylated amine (N connected to C(=O)C)
    acetyl_amide = Chem.MolFromSmarts('[NH]C(=O)C')
    if not mol.HasSubstructMatch(acetyl_amide):
        return False, "No acetylated amine group"
    
    return True, "Contains carboxylic acid and acetylated amine group"