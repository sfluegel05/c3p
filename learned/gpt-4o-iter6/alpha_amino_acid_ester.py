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

    # Relaxed pattern for alpha-amino acid to allow for more structures
    # Example: find NR-C(R)(R)-C(=O)O where N is connected to central C and there's an ester
    alpha_amino_acid_ester_pattern = Chem.MolFromSmarts("[NX3][CX4](R)[C](=O)O[!H]")
    
    if mol.HasSubstructMatch(alpha_amino_acid_ester_pattern):
        return True, "Contains likely alpha-amino acid backbone with ester linkage"
    
    return False, "Does not match common alpha-amino acid ester pattern"