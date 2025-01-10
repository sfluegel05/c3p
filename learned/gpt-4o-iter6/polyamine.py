"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic compound with two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to match amino groups
    # Matches primary (NH2), secondary (NHR), tertiary (NR2) amines
    # Ensure we exclude non-amino nitrogen groups such as nitro groups
    primary_amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O),!$(N~N)]")  # NH2-R or NH-R
    secondary_tertiary_amino_pattern = Chem.MolFromSmarts("[NX3;H0;!$(NC=O),!$(N~N)]")  # NR2-R (tri/quaternary)
    
    # Find all matches for amino groups in the molecule
    primary_amino_matches = mol.GetSubstructMatches(primary_amino_pattern)
    secondary_tertiary_amino_matches = mol.GetSubstructMatches(secondary_tertiary_amino_pattern)
    
    # Combine both matches to find unique amino groups
    unique_amino_matches = set(primary_amino_matches) | set(secondary_tertiary_amino_matches)
    num_amino_groups = len(unique_amino_matches)

    # Check if there are 2 or more amino groups
    if num_amino_groups >= 2:
        return True, f"Contains {num_amino_groups} amino groups"
    else:
        return False, f"Only {num_amino_groups} amino group(s) found, need at least 2"