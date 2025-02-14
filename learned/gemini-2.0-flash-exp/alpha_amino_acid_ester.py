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

    # SMARTS pattern for alpha carbon with amino group (allowing for any other attached groups)
    alpha_amino_pattern = Chem.MolFromSmarts("[CX4H0-2]([NX3H0-2])([CX4,CX3])") 
    # SMARTS pattern for any ester group
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
    #SMARTS pattern for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")

    #check that the patterns are matched
    alpha_amino_matches = mol.GetSubstructMatches(alpha_amino_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not alpha_amino_matches:
        return False, "Molecule does not contain an alpha-amino acid core"
    
    if not ester_matches:
         return False, "Molecule does not contain an ester group"

    
    if carboxylic_acid_matches:
        return False, "Molecule contains a free carboxylic acid, not an ester"

    return True, "Molecule contains an alpha-amino acid with an ester group."