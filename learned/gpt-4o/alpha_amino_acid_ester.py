"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is characterized by an amine (NH2) group at the alpha carbon
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

    # Patterns for alpha-amino acid core:
    # Capture amine group on alpha carbon (N-C-C(=O)) where C(=O) is esterified
    amino_acid_pattern = Chem.MolFromSmarts("[$([NX3;H2,H1;!$(NC=O)])]-[C;H2,H1]([C]=O)-[O]!@[CX3](=O)")

    # Patterns for ester linkage - Esterification should be directly on the carboxylic group
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[O;!H1][C]")

    # Verify presence of an alpha-amino structure
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha-amino acid substructure found"
    
    # Verify presence of ester linkage
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    return True, "Contains an alpha-amino acid substructure with ester linkages"