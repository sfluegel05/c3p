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
    
    # General pattern for L-alpha-amino acid backbone with stereochemistry
    l_alpha_amino_acid_pattern = Chem.MolFromSmarts('[NX3;H2][C@H]([CX4])[CX3](=O)O')  # Primary amine with chiral alpha-carbon and carboxyl
    if not mol.HasSubstructMatch(l_alpha_amino_acid_pattern):
        return (False, "No L-alpha-amino acid backbone found")

    # Pattern for nitrogen with a connected acyl group
    n_acyl_pattern = Chem.MolFromSmarts('[NX3;H1][CX3](=O)[!#1]')
    n_acyl_matches = mol.GetSubstructMatches(n_acyl_pattern)

    if not n_acyl_matches:
        return (False, "No N-acyl group found attached to nitrogen")

    # Locate if the N-acyl is attached to the nitrogen within the context of an L-alpha-amino acid pattern
    for match in n_acyl_matches:
        nitrogen_idx = match[0]
        for neighbor in mol.GetAtomWithIdx(nitrogen_idx).GetNeighbors():
            idx = neighbor.GetIdx()
            # The alpha carbon should be part of the L-alpha-amino acid structure and connected to the n-acyl nitrogen
            if mol.HasSubstructMatch(l_alpha_amino_acid_pattern, atoms=[idx]):
                return (True, "Contains L-alpha-amino acid backbone with an N-acyl substituent correctly attached")

    return (False, "N-acyl group not correctly integrated into L-alpha-amino acid nitrogen")