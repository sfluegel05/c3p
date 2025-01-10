"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group attached to the nitrogen of an amino acid backbone.

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

    # Check for N-acetyl group attached to nitrogen
    n_acetyl_pattern = Chem.MolFromSmarts("N(C(=O)C)")
    acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    if not acetyl_matches:
        return False, "No N-acetyl group attached to nitrogen atom"

    # Check for amino acid backbone
    # Alpha carbon connected to nitrogen and carboxyl group
    amino_acid_pattern = Chem.MolFromSmarts("[CH1](N(C(=O)C))C(=O)[O,H]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "No amino acid backbone found"

    # Optional: Ensure the molecule is an alpha-amino acid (one central carbon)
    # and the acetyl group is specifically attached to the nitrogen of the amino group
    # Additional checks can be implemented if necessary

    return True, "Molecule is an N-acetyl-amino acid"