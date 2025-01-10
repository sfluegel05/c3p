"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify amino groups (primary, secondary, and possibly tertiary amines)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Count the amino groups
    num_amino_groups = len(amino_matches)
    
    if num_amino_groups >= 2:
        return True, f"Molecule contains {num_amino_groups} amino groups, indicating a polyamine"
    else:
        return False, f"Molecule contains {num_amino_groups} amino group(s), fewer than required for a polyamine"