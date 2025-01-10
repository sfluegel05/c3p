"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid contains an acetyl group bound to an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for Acetyl group pattern "CC(=O)N"
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group found"
    
    # Look for amino acid moiety pattern "NC(C(O)=O)"
    amino_acid_pattern = Chem.MolFromSmarts("NC(C(O)=O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid moiety found"
    
    # Ensure acetyl group is connected to the amino acid moiety
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    for acetyl in acetyl_matches:
        for amino_acid in amino_acid_matches:
            # Check if the nitrogen in acetyl is the same as the nitrogen in amino acid
            # This checks connectivity
            if acetyl[2] == amino_acid[0]:
                return True, "Contains acetyl group bound to amino acid"
    
    return False, "Acetyl group not directly connected to amino acid moiety"