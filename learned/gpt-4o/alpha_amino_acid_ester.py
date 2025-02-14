"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is derived from the formal condensation of an alpha-amino acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define alpha-amino acid ester pattern:
    # N-C-C(=O)O-R (where CH can be chiral)
    amino_acid_ester_pattern = Chem.MolFromSmarts("N[C;!R][C](=O)O")
    
    # Check for presence of the pattern
    if mol.HasSubstructMatch(amino_acid_ester_pattern):
        return True, "Contains the alpha-amino acid ester functional group"
    
    return False, "Does not contain the alpha-amino acid ester functional group"