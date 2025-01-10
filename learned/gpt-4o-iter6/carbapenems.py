"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determine if a given SMILES string is a carbapenem.

    Args:
        smiles (str): The SMILES string of the chemical entity.
    
    Returns:
        bool: True if SMILES represents a carbapenem, False otherwise.
        str: Explanation of the reasoning for the classification.
    """
    
    # Define the refined SMARTS pattern for the carbapenem bicyclic core:
    # Structure with a beta-lactam ring and common substitution positions
    core_pattern = Chem.MolFromSmarts("O=C1N2C([C@@H]3C[C@@H]23)C1")  # Bicyclic [3.2.0] with a beta-lactam ring

    # Additional patterns capturing functional groups and typical substitutions
    additional_patterns = [
        Chem.MolFromSmarts("C(=O)O"),  # Carboxylate group, usually present
        Chem.MolFromSmarts("CS"),  # Thioether group, common in side chains
        Chem.MolFromSmarts("CC(=O)N"),  # Amide linkage
    ]

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for the core carbapenem pattern in the molecule
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match the carbapenem core structure"

    # Check for additional common substructures
    missing_structures = []
    for pattern in additional_patterns:
        if not mol.HasSubstructMatch(pattern):
            missing_structures.append("Missing functional group matching pattern: " + Chem.MolToSmiles(pattern))
    if missing_structures:
        return False, "; ".join(missing_structures)

    # If core and additional features are present, classify as carbapenem
    return True, "Contains characteristic carbapenem core and common functional groups"