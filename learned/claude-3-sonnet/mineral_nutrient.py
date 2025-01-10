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
        'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn',  # Major minerals
        'Ba', 'Sr', 'La', 'Cs', 'Al', 'Sb'   # Less common but nutritionally relevant
    }
    
    # Common nutritionally relevant anions and their components
    common_anions = {
        'Cl', 'F', 'P', 'S', 'Si', 'N',  # Inorganic
        'C', 'O', 'H'  # For carbonates, phosphates, hydroxides etc.
    }

    # Get all elements and charges present
    elements = set()
    has_cation = False
    has_anion = False
    for atom in mol.GetAtoms():
        elements.add(atom.GetSymbol())
        charge = atom.GetFormalCharge()
        if charge > 0:
            has_cation = True
        elif charge < 0:
            has_anion = True

    # Check if at least one essential mineral cation is present
    essential_cations_found = elements.intersection(essential_cations)
    if not essential_cations_found:
        return False, "No essential mineral cations found"

    # Handle special cases of non-ionic mineral compounds
    direct_metal_bonds = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in essential_cations:
            if len(atom.GetBonds()) > 0:
                direct_metal_bonds = True

    # Count atoms to identify molecular complexity
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H')
    
    # Verify ionic nature or acceptable direct metal bonds
    if not (has_cation or has_anion) and not direct_metal_bonds:
        return False, "Neither ionic nor valid metal compound"

    # Count water molecules
    water_pattern = Chem.MolFromSmarts("[H]O[H]")
    if water_pattern:
        water_matches = len(mol.GetSubstructMatches(water_pattern))
    else:
        water_matches = 0
        
    # Exclude complex organic compounds
    if carbon_count > 2:  # Allow only simple carbonates, bicarbonates
        if not (carbon_count == oxygen_count and has_cation):  # Exception for carbonates
            return False, "Too complex - not a simple mineral nutrient"

    # Check for rings (excluding special cases)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        return False, "Contains rings - too complex for mineral nutrient"

    # Check remaining elements are acceptable anions
    non_essential_elements = elements - essential_cations
    if not non_essential_elements.issubset(common_anions):
        unexpected_elements = non_essential_elements - common_anions
        return False, f"Contains unexpected elements: {unexpected_elements}"

    # Build reason string
    reason_parts = []
    if direct_metal_bonds:
        reason_parts.append("metal-containing compound")
    elif has_cation or has_anion:
        reason_parts.append("ionic mineral compound")
    
    reason_parts.append(f"contains essential elements: {', '.join(essential_cations_found)}")
    
    if water_matches > 0:
        reason_parts.append(f"hydrated form ({water_matches} water molecules)")
    
    return True, " - ".join(reason_parts)