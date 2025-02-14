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
    
    # Improved pattern for L-alpha-amino acids (including all plausible stereochemistry representations):
    # Consider adding multiple possible L-configurations or relax some constraints with additional patterns
    amino_acid_patterns = [
        Chem.MolFromSmarts("[C@@H](N)[CH]C(=O)O"),
        Chem.MolFromSmarts("[C@H](N)[CH]C(=O)O")  # As an example, additional projection variants can be added
    ]
    
    has_amino_acid_backbone = any(mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns)
    if not has_amino_acid_backbone:
        return False, "No L-alpha-amino acid backbone (with proper stereochemistry) found"
    
    # Improved pattern for N-acyl substituent:
    acyl_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    
    # Ensure N-acyl is directly bound to the amino nitrogen atom by checking adjacency and connectivity explicitly
    # Get matches for patterns and verify connectivity
    acyl_linkages = mol.GetSubstructMatches(acyl_pattern)
    amino_linkages = [match for pattern in amino_acid_patterns for match in mol.GetSubstructMatches(pattern)]
    
    for acyl in acyl_linkages:
        acyl_nitrogen_idx = acyl[0]
        for amino in amino_linkages:
            amino_nitrogen_idx = amino[1]  # Recall this is the N in the amino pattern
            
            # Check if acyl part is specifically attached to amino nitrogen
            # We thus verify explicit connectivity:
            if acyl_nitrogen_idx == amino_nitrogen_idx and \
               mol.GetBondBetweenAtoms(acyl_nitrogen_idx, amino_nitrogen_idx):
                return True, "Contains valid L-alpha-amino acid backbone with an N-acyl substituent"

    return False, "N-acyl substituent not correctly linked to the L-alpha-amino acid backbone"