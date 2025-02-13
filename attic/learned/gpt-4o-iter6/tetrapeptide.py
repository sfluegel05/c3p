"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide must contain four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized peptide bond pattern N-C(=O)
    peptide_pattern = Chem.MolFromSmarts("[NX3,NX2H1,NX2H2][C](=O)")  # captures peptide/cyclic bonds
    
    # Checking at least 3 peptide linkages implies a tetrapeptide
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) < 3:
        return False, f"Contains {len(peptide_matches)} peptide bonds, expected at least 3"

    # Refine the amino acid backbone pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX2H1,NX2H2][CX4][CX3](=O)")  # generalized backbone with flexibility
    
    # Look for the presence of four amino residues
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) < 4:
        return False, f"Detected {len(amino_acid_matches)} amino acid residues, expected 4"

    return True, "Contains four amino-acid residues connected by peptide linkages"