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
    
    # Refined SMARTS pattern to match amino groups
    # Includes primary (NH2), secondary (NHR), tertiary (NR2) amines,
    # and quaternary ammonium groups
    amino_pattern = Chem.MolFromSmarts("[NX3,NX4+;!$(NC=O),!$(N~N)]")

    # Find all matches for amino groups in the molecule
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Check if there are 2 or more amino groups
    num_amino_groups = len(amino_matches)
    if num_amino_groups >= 2:
        return True, f"Contains {num_amino_groups} amino groups"
    else:
        return False, f"Only {num_amino_groups} amino group(s) found, need at least 2"