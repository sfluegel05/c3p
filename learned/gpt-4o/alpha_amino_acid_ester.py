"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is characterized by an amino group (NH2) at the alpha carbon
    of the carboxylic acid, which is esterified to an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update pattern to capture an alpha-amino acid core more reliably
    # This pattern captures NH2 on alpha-carbon directly bonded to a carboxylic group
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]-[C;H2,H1](-[C;H1]C(=O)[OX2H0]-[CX4])[CX3](=O)")

    # Check if molecule matches the alpha-amino acid ester structure
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha-amino acid ester substructure found"
    
    return True, "Contains an alpha-amino acid core with ester linkages"