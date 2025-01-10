"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a compound is a mineral nutrient based on its SMILES string.
    
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

    # Essential metallic elements typically found in nutrients
    essential_metallic_elements = {'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn', 'Ba', 'La', 'Cs', 'Al', 'Sb', 'Pd'}
    # Nutrient anions indicative of common inorganic nutrient minerals
    nutrient_anions_smarts = [
        "[Cl]", "[F]", "[O-]P([O-])(=O)", "[O-]C([O-])=O",  # Typical inorganic anions
        "[O-]S([O-])(=O)"]  # Sulfate or Phosphate derivatives
    ]

    # Check for the presence of at least one essential metallic element
    metal_found = any(atom.GetSymbol() in essential_metallic_elements for atom in mol.GetAtoms())
    if not metal_found:
        return False, "No essential metallic elements found"

    # Ensure relevant inorganic anions for nutrient character are present
    anion_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(anion)) for anion in nutrient_anions_smarts)
    if not anion_found:
        return False, "No relevant mineral nutrient anions found"

    # Allow for longer carbon chains if balanced by the presence of magnesium, calcium, etc.
    # Count carbon atoms but condone if present with significant nutrient metallic atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count > 36:
        return False, f"Contains {carbon_count} carbon atoms suggesting an overly complex organic molecule"

    return True, "Structure consistent with mineral nutrients containing essential elements and simple inorganic anions"