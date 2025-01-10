"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:26144 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic compound containing two or more amino groups.

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

    # Check if molecule is organic (contains carbon)
    has_carbon = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            has_carbon = True
            break
    if not has_carbon:
        return False, "Not an organic compound"

    # Count different types of amino groups
    
    # Primary amines (-NH2)
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2]")
    primary_amines = len(mol.GetSubstructMatches(primary_amine_pattern))
    
    # Secondary amines (-NH-)
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H1]")
    secondary_amines = len(mol.GetSubstructMatches(secondary_amine_pattern))
    
    # Tertiary amines (-N<)
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3H0]")
    tertiary_amines = len(mol.GetSubstructMatches(tertiary_amine_pattern))
    
    # Protonated amines (NH3+, NH2+)
    protonated_amine_pattern = Chem.MolFromSmarts("[NX4+]")
    protonated_amines = len(mol.GetSubstructMatches(protonated_amine_pattern))

    # Count total amine groups
    total_amines = primary_amines + secondary_amines + tertiary_amines + protonated_amines

    # Exclude certain nitrogen-containing groups that shouldn't count
    
    # Nitro groups
    nitro_pattern = Chem.MolFromSmarts("[$([NX3](=O)=O),$([NX3+](=O)[O-])]")
    nitro_groups = len(mol.GetSubstructMatches(nitro_pattern))
    
    # Amide groups
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_groups = len(mol.GetSubstructMatches(amide_pattern))
    
    # Subtract these from total amine count
    true_amine_count = total_amines - nitro_groups - amide_groups

    if true_amine_count < 2:
        return False, f"Contains only {true_amine_count} amino groups, minimum 2 required"

    # Success case
    return True, f"Contains {true_amine_count} amino groups (Primary: {primary_amines}, Secondary: {secondary_amines}, Tertiary: {tertiary_amines}, Protonated: {protonated_amines})"