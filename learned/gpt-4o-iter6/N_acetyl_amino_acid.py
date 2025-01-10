"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is defined as having an acetyl group directly connected to an amino nitrogen
    which is part of an amino acid structure.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an N-acetyl-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acetyl group attached to nitrogen in form: CC(=O)N
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    if not acetyl_matches:
        return False, "No acetyl group directly bonded to nitrogen found"

    # General pattern for an amino acid: N[C;!R][CX4]C(=O)O (allow flexibility in R-group)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4,CX3][CX3](=O)[OX1-,OX2H1]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "No amino acid moiety found"
    
    # Ensure acetyl and amino acid patterns are part of the same structure
    for acetyl in acetyl_matches:
        acetyl_nitrogen = acetyl[2]  # Nitrogen from acetyl group
        for amino_acid in amino_acid_matches:
            amino_acid_nitrogen = amino_acid[0]  
            
            # Check if the acetyl nitrogen is the same as the amino acid nitrogen
            if acetyl_nitrogen == amino_acid_nitrogen:
                return True, "Contains acetyl group directly connected to an amino acid moiety"
    
    return False, "Acetyl group not directly connected to the correct amino acid moiety"