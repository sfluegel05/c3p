"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate a list of common amino acid SMILES patterns to identify amino acids in a polypeptide
    # These are simplified, representative substructures of amino acids at the amine or carboxyl terminus.
    amino_acid_patterns = [
        Chem.MolFromSmarts('[NX3][CX4](C)C=O'),  # Glycine pattern, extends to other amino acids with similar backbone
        Chem.MolFromSmarts('[NX3][CX4H2]CC=O'),  # Alanine pattern
        Chem.MolFromSmarts('[NX3][CX4](CC)[CX3](=O)'),  # Random extended pattern to cover potential modifications
        # Add more patterns or variations for common and modified amino acids
    ]

    # Function to match and count amino acid patterns in the given molecule
    num_amino_acids = 0
    for pattern in amino_acid_patterns:
        matches = mol.GetSubstructMatches(pattern)
        num_amino_acids += len(matches)

    # Consider each peptide bond (represented correctly) represents two amino acids, adjust as necessary
    # This is a simplification, real molecules may need closer inspection
    # Typically each peptide bond implies the presences of two amino acids in chain, so half the simple counts
    if num_amino_acids / 2 >= 10:
        return True, f"Contains approximately {num_amino_acids // 2} amino acid residues, classifying as polypeptide"
    else:
        return False, f"Contains approximately {num_amino_acids // 2} amino acid residues, not enough for a polypeptide"