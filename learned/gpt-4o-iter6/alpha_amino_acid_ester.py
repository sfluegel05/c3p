"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester should have an alpha-carbon (with both an amino and a carboxyl group)
    and a carboxyl group that is esterified.

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

    # Define the alpha-amino acid ester pattern
    # -N[C@H](C)C(=O)O linking to an alkyl or aryl group for the ester

    # Find carboxyl linked to ester (-C(=O)O)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Find alpha-amino acid pattern with ester linkage
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)][C@H](C)C(=O)O")
    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Contains alpha-amino acid backbone with ester linkage"
    
    return False, "Does not match alpha-amino acid ester pattern"