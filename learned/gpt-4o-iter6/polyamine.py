"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic compound with two or more amino groups, including different hydrogenation or charged states.

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
    
    # Refined SMARTS pattern to match amino groups, including both charged and neutral amines
    # Matches primary (NH2), secondary (NH), tertiary (NR) amines, quaternary amines and their charged counterparts
    amine_pattern = Chem.MolFromSmarts("[NX3,NX4+;!$(NC=O)]")  # Ensure non-amide nitrogen

    # Find matches for amino groups in the molecule
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    num_amino_groups = len(amine_matches)

    # Check if there are 2 or more amino groups
    if num_amino_groups >= 2:
        return True, f"Contains {num_amino_groups} amino groups"
    else:
        return False, f"Only {num_amino_groups} amino group(s) found, need at least 2"