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
    
    # Identify primary amine pattern (should capture NH2 groups)
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2]")
    primary_amine_matches = mol.GetSubstructMatches(primary_amine_pattern)

    # Identify secondary amine pattern (focus on NH with two carbons attached)
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H1]([CX4])[CX4]")
    secondary_amine_matches = mol.GetSubstructMatches(secondary_amine_pattern)

    # Identify charged nitrogen which might indicate amino groups in ionic form
    charged_amine_pattern = Chem.MolFromSmarts("[NX4+]")
    charged_amine_matches = mol.GetSubstructMatches(charged_amine_pattern)
    
    # Total possible amino group matches 
    num_amino_groups = len(primary_amine_matches) + len(secondary_amine_matches) + len(charged_amine_matches)
    
    if num_amino_groups >= 2:
        return True, f"Molecule contains {num_amino_groups} amino groups, indicating a polyamine"
    else:
        return False, f"Molecule contains {num_amino_groups} amino group(s), fewer than required for a polyamine"