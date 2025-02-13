"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    Requires both a stereospecific L-alpha-amino acid backbone and an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved pattern for L-alpha-amino acids (include stereophysical details ensuring L-config):
    amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)[CH]C(=O)O")  # General pattern with focus on L-stereochemistry
    
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone (with proper stereochemistry) found"

    # Improved pattern for N-acyl substituent:
    acyl_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")

    # Ensure N-acyl is directly bound to the amino nitrogen atom
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No N-acyl substituent found attached to the Nitrogen of the amino acid"

    # Collecting additional evidence on exact positioning
    # Assume it must check if N-acyl is a substituent of the nitrogen in the core alpha-amino acid
    acyl_linkages = mol.GetSubstructMatches(acyl_pattern)
    amino_linkages = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Check if the acyl part is specifically attached to the L-alpha amino component of nitrogen
    for acyl in acyl_linkages:
        for amino in amino_linkages:
            # Ensure the N of N-acyl is the same as the N in the amino acid pattern
            if mol.GetAtomWithIdx(acyl[0]).GetIdx() == mol.GetAtomWithIdx(amino[1]).GetIdx():
                return True, "Contains valid L-alpha-amino acid backbone with an N-acyl substituent"

    return False, "N-acyl substituent not correctly linked to the L-alpha-amino acid backbone"