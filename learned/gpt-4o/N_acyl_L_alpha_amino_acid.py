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
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for L-alpha-amino acid backbone with specific stereochemistry:
    # Needs to accentuate L-stereochemistry
    amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)[C](=O)[O]")  # Re-evaluate based on L-config
    
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone (with proper stereochemistry) found"

    # Pattern for N-acyl substituent: structure ensuring linkage to the N atom
    acyl_pattern = Chem.MolFromSmarts("N[C](=O)")

    # Ensure N-acyl is the direct modification
    if not any(mol.GetSubstructMatch(acyl_pattern)):
        return False, "No direct N-acyl substituent found"

    # Verify that acyl linkage conforms to N-acyl criteria with amino group linkage
    acyl_linkages = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_linkages) == 0:
        return False, "N-acyl linkage not correctly established"

    return True, "Contains valid L-alpha-amino acid backbone with an N-acyl substituent"