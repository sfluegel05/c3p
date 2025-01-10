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
        return None, "Invalid SMILES string"
    
    # SMARTS patterns for primary and secondary amines, and charged amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # captures primary amines not bound to carbonyl
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H][!,$([#6]=O)]")  # captures NH attached to two C (not amidic)
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([C])[C]")  # captures three Cs attached
    charged_amine_pattern = Chem.MolFromSmarts("[NX4+]")
    
    # Identifying all amino group matches
    primary_amine_matches = mol.GetSubstructMatches(primary_amine_pattern)
    secondary_amine_matches = mol.GetSubstructMatches(secondary_amine_pattern)
    tertiary_amine_matches = mol.GetSubstructMatches(tertiary_amine_pattern)
    charged_amine_matches = mol.GetSubstructMatches(charged_amine_pattern)
    
    # Aggregate total number of amine groups found
    num_amino_groups = (len(primary_amine_matches) + 
                        len(secondary_amine_matches) + 
                        len(tertiary_amine_matches) +
                        len(charged_amine_matches))
    
    if num_amino_groups >= 2:
        return True, f"Molecule contains {num_amino_groups} amino groups, indicating a polyamine"
    else:
        return False, f"Molecule contains {num_amino_groups} amino group(s), fewer than required for a polyamine"