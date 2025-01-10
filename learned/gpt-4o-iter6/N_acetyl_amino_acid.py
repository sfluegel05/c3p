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
    
    # Look for acetyl group attached to nitrogen
    acetyl_amino_pattern = Chem.MolFromSmarts("CC(=O)N[*]")
    if not mol.HasSubstructMatch(acetyl_amino_pattern):
        return False, "No acetyl group directly bonded to the nitrogen of an amino group"
    
    # General pattern for amino acids
    # This includes flexibility in backbone connections
    amino_acid_backbone_pattern = Chem.MolFromSmarts("N[C;!R]C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_backbone_pattern):
        return False, "No amino acid moiety found"
    
    # Ensure acetyl and amino acid patterns are part of the same structure
    acetyl_matches = mol.GetSubstructMatches(acetyl_amino_pattern)
    backbone_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)
    for acetyl in acetyl_matches:
        for backbone in backbone_matches:
            # Check if the acetyl nitrogen is the same as the amino acid nitrogen
            if acetyl[2] == backbone[0]:
                return True, "Contains acetyl group directly connected to an amino acid moiety"
    
    return False, "Acetyl group not directly connected to the correct amino acid moiety"