"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a compound is a mineral nutrient based on its SMILES string.
    Focuses on essential inorganic nutrients with simple structures.
    
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

    # Define essential metallic elements typically found in nutrients
    essential_metallic_elements = {'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn'}
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

    # Ensure the structure is relatively simple (e.g., not too many rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Complex organic structures detected"

    # Determine the simplicity of structure (controlled number of elements)
    if mol.GetNumAtoms() > 25:
        return False, "Structure too complex to be a typical nutrient"

    return True, "Structure consistent with mineral nutrients containing essential elements and simple inorganic anions"