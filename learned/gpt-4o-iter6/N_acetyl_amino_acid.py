"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is defined as having an acetyl group directly connected to an amino nitrogen
    which is part of an amino acid structure. This function allows flexibility to cover various amino acids.
    
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
    acetyl_to_nitrogen_pattern = Chem.MolFromSmarts("CC(=O)N")
    acetyl_to_nitrogen_matches = mol.GetSubstructMatches(acetyl_to_nitrogen_pattern)
    if not acetyl_to_nitrogen_matches:
        return False, "No acetyl group directly bonded to nitrogen found"

    # General pattern for an amino acid: N-R-C-C(=O)O, allowing diverse R-group variants.
    # Accommodate for different variations like L/D stereochemistry and different R-groups.
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[OX1-,OX2H1]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "No amino acid moiety found"
    
    # Ensure acetyl nitrogen is the part of amino acid nitrogen
    for acetyl in acetyl_to_nitrogen_matches:
        acetyl_nitrogen = acetyl[2]  # Nitrogen from the acetyl group
        for amino_acid in amino_acid_matches:
            # check if acetyl nitrogen is part of the amino acid structure's nitrogen
            if acetyl_nitrogen == amino_acid[0]:
                return True, "Contains acetyl group directly connected to an amino acid moiety"
    
    return False, "Acetyl group not directly connected to the correct amino acid moiety"