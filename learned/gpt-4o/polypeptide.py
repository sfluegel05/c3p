"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved peptide bond pattern: N-C(=O)-C(alpha carbon of amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C][C](=O)C")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Count the number of matched peptide motifs
    num_amino_acid_residues = len(peptide_matches)

    # Polypeptide criteria: 10 or more residues
    if num_amino_acid_residues >= 10:
        return True, f"Contains {num_amino_acid_residues} amino acid residues, qualifies as a polypeptide"
    else:
        return False, f"Contains only {num_amino_acid_residues} amino acid residues, not a polypeptide"

    return None, None