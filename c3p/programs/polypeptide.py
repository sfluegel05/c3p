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
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined pattern to match a broader range of amino acid structures
    # This pattern captures diverse amino acids, including standard and some modified residues
    amino_acid_patterns = [
        Chem.MolFromSmarts('[NX3][C@@H](C)C(=O)O'),  # Standard amino acid pattern
        Chem.MolFromSmarts('[NX3][C@@H]([C@@H](C)O)C(=O)O'),  # Example for hydroxy-amino acids, etc.
        # Further patterns can be added to capture other common modifications
    ]

    num_amino_acids = 0
    
    for pattern in amino_acid_patterns:
        matches = mol.GetSubstructMatches(pattern)
        num_amino_acids += len(matches)

    if num_amino_acids >= 10:
        return True, f"Contains {num_amino_acids} amino acid residues, classifying as polypeptide"
    else:
        return False, f"Contains {num_amino_acids} amino acid residues, not enough for a polypeptide"