"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    Mineral nutrients are inorganic substances required for essential metabolic functions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Essential mineral nutrient cations
    essential_cations = {
        'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn', 'Cu', 'Mn',  # Major and trace minerals
        'Ba', 'Sr', 'La', 'Cs', 'Al',  # Less common but nutritionally relevant
        'Pd', 'Sb'  # Special cases
    }
    
    # Common nutritionally relevant anions
    common_anions = {
        'Cl', 'F', 'P', 'S', 'Si',  # Inorganic
        'C', 'O'  # For carbonates, phosphates, etc.
    }

    # Get all elements and charges present
    elements = set()
    has_cation = False
    has_anion = False
    for atom in mol.GetAtoms():
        elements.add(atom.GetSymbol())
        if atom.GetFormalCharge() > 0:
            has_cation = True
        elif atom.GetFormalCharge() < 0:
            has_anion = True

    # Check if at least one essential mineral cation is present
    essential_cations_found = elements.intersection(essential_cations)
    if not essential_cations_found:
        return False, "No essential mineral cations found"

    # Verify ionic nature for most compounds
    if not (has_cation or has_anion) and not any(x in elements for x in ['La', 'Sb']):
        return False, "Not an ionic compound"

    # Count atoms to identify molecular complexity
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    # Allow higher carbon counts only for specific nutritional organic salts
    if carbon_count > 20:  # Allow for stearates and similar long-chain salts
        if not (has_cation and has_anion and oxygen_count >= 2):
            return False, "Organic compound too complex for mineral nutrient"
    
    # Check for complex structures
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        return False, "Contains rings - too complex for mineral nutrient"

    # Check remaining elements are acceptable anions
    non_essential_elements = elements - essential_cations
    if not non_essential_elements.issubset(common_anions):
        unexpected_elements = non_essential_elements - common_anions
        return False, f"Contains unexpected elements: {unexpected_elements}"

    # Specific exclusions for toxic or non-nutritive forms
    if 'Se' in elements and carbon_count > 0:  # Exclude organic selenium compounds
        return False, "Organic selenium compounds not considered mineral nutrients"
    
    if 'S' in elements and not has_cation:  # Exclude pure sulfur compounds
        return False, "Pure sulfur compounds not considered mineral nutrients"

    # Success case - construct reason string
    reason_parts = []
    if has_cation or has_anion:
        reason_parts.append("ionic compound")
    reason_parts.append(f"contains essential elements: {', '.join(essential_cations_found)}")
    
    if carbon_count > 3:
        reason_parts.append("nutritionally relevant organic salt")
    elif carbon_count > 0:
        reason_parts.append("simple inorganic/organic salt")
    else:
        reason_parts.append("inorganic nutrient")
    
    return True, " - ".join(reason_parts)