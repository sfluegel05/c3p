"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a compound is a mineral nutrient based on its SMILES string.
    Mineral nutrients are typically inorganic compounds containing metallic elements
    necessary for metabolic/structural functions in the human body.
    
    Args:
        smiles (str): SMILES string of the compound
        
    Returns:
        bool: True if compound is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define common metallic and essential non-metallic elements in mineral nutrients
    metallic_elements = {'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn', 'Al', 'Ba', 'Cs', 'Sb', 'Pd', 'La'}
    non_metallic_anions = {'P', 'S', 'Cl', 'F', 'O', 'N'}  # Part of common inorganic anions

    # Check for presence of at least one metallic element
    metal_found = any(Chem.MolFromSmarts(f"[#{elem}]") and mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{elem}]")) for elem in metallic_elements)
    if not metal_found:
        return False, "No metallic elements found"

    # Check for presence of common inorganic anions (e.g., sulfates, phosphates, chlorides)
    anion_found = any(atom.GetSymbol() in non_metallic_anions for atom in mol.GetAtoms())
    if not anion_found:
        return False, "No common inorganic anions or ions found"

    return True, "Contains metallic elements and common inorganic anions indicative of a mineral nutrient"