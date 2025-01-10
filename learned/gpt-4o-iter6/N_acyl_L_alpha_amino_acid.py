"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    This class is defined by an L-alpha-amino acid backbone with an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # General pattern for L-alpha-amino acid backbone
    amino_acid_pattern = Chem.MolFromSmarts('[NX3][C@H]([*])[C](=O)O')  # Chiral carbon attached to nitrogen and carboxyl group
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return (False, "No L-alpha-amino acid backbone found")

    # Pattern for N-acyl group connected to the alpha nitrogen atom
    n_acyl_pattern = Chem.MolFromSmarts('[NX3][C](=O)[*]')  # N-acyl connected to nitrogen atom
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return (False, "No N-acyl group found")

    # Additionally, ensure the N-acyl group is directly attached to the amino acid nitrogen
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    for match in matches:
        # Confirm the acyl group is connected to the same nitrogen
        nitrogen_idx = next(idx for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == 'N')
        n_acyl_matches = [m for m in mol.GetSubstructMatches(n_acyl_pattern) if nitrogen_idx in m]
        if n_acyl_matches:
            return (True, "Contains L-alpha-amino acid backbone with an N-acyl substituent")

    return (False, "N-acyl group not correctly attached to L-alpha-amino acid nitrogen")