"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic compound with two or more amino groups, including their charged counterparts.

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
    
    # SMARTS pattern to match amino groups, including charged states
    # Matches primary (NH2), secondary (NH), tertiary (NR) amines and charged amines
    primary_amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;+0]")
    secondary_amino_pattern = Chem.MolFromSmarts("[NX3;H1;+0]")
    tertiary_amino_pattern = Chem.MolFromSmarts("[NX3;H0;+0]")
    charged_amino_pattern = Chem.MolFromSmarts("[NX3;+1]")
    
    # Find all matches for amino groups in the molecule
    primary_amino_matches = mol.GetSubstructMatches(primary_amino_pattern)
    secondary_amino_matches = mol.GetSubstructMatches(secondary_amino_pattern)
    tertiary_amino_matches = mol.GetSubstructMatches(tertiary_amino_pattern)
    charged_amino_matches = mol.GetSubstructMatches(charged_amino_pattern)
    
    # Combine all matches to find unique amino groups, considering their charge
    unique_amino_matches = set(primary_amino_matches) | set(secondary_amino_matches) | set(tertiary_amino_matches)
    unique_amino_matches.update(charged_amino_matches)
    num_amino_groups = len(unique_amino_matches)

    # Check if there are 2 or more amino groups
    if num_amino_groups >= 2:
        return True, f"Contains {num_amino_groups} amino groups"
    else:
        return False, f"Only {num_amino_groups} amino group(s) found, need at least 2"