"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a compound is a mineral nutrient based on its SMILES string.
    Extends the range of essential inorganic nutrients with a broader set of elements.
    
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

    # Define broader range of essential metallic elements typically found in nutrients
    essential_metallic_elements = {'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn', 'Ba', 'La', 'Cs', 'Al', 'Sb', 'Pd'}
    # Define essential non-metallic elements associated with nutrients
    nutrient_anions = {'P', 'S', 'Cl', 'F', 'O', 'N'}  # Common in simple inorganic salts

    # Check for presence of at least one essential metallic element
    metal_found = any(atom.GetSymbol() in essential_metallic_elements for atom in mol.GetAtoms())
    if not metal_found:
        return False, "No essential metallic elements found"

    # Check for presence of common inorganic anions (e.g., sulfates, phosphates, chlorides) in simpler structures
    anion_found = any(atom.GetSymbol() in nutrient_anions for atom in mol.GetAtoms())
    if not anion_found:
        return False, "No relevant mineral nutrient anions found"

    # Allow more complex structures since some valid nutrients might include them
    # for magnesium distearate or other larger salts/hydrates
    if mol.GetNumAtoms() > 50:
        return False, "Structure unnecessarily complex to be a typical nutrient"

    return True, "Structure consistent with mineral nutrients containing essential elements and simple inorganic anions"