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

    # Essential mineral elements (including trace elements)
    essential_elements = {
        'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn', 'Cu', 'Mn', 'I', 'Se', 'Mo',
        'Cr', 'Co', 'F', 'P', 'Cl', 'S', 'Al', 'Ba', 'Sr', 'La', 'Cs',
        'Sb', 'Pd'  # Including some less common but potentially beneficial elements
    }
    
    # Common anions in mineral nutrients
    common_anions = {'Cl', 'F', 'P', 'S', 'O', 'C', 'N', 'Si'}

    # Get all elements present in the molecule
    elements = set()
    for atom in mol.GetAtoms():
        elements.add(atom.GetSymbol())

    # Check if at least one essential mineral element is present
    essential_elements_found = elements.intersection(essential_elements)
    if not essential_elements_found:
        return False, "No essential mineral elements found"

    # Count carbon atoms - too many carbons suggest organic compound
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count > 3:  # Allow simple carbonates, acetates, etc.
        return False, "Too many carbon atoms for a mineral nutrient"

    # Check for complex organic structures
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        return False, "Contains rings - too complex for mineral nutrient"

    # Look for ionic species or simple covalent compounds
    has_ionic_species = False
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            has_ionic_species = True
            break
    
    # Check if remaining elements are common anions
    non_essential_elements = elements - essential_elements
    if not non_essential_elements.issubset(common_anions):
        unexpected_elements = non_essential_elements - common_anions
        return False, f"Contains unexpected elements: {unexpected_elements}"

    # Additional checks for molecular complexity
    if rdMolDescriptors.CalcNumRotatableBonds(mol) > 5:
        return False, "Structure too complex for mineral nutrient"

    # Success case - construct reason string
    reason_parts = []
    if has_ionic_species:
        reason_parts.append("ionic compound")
    reason_parts.append(f"contains essential elements: {', '.join(essential_elements_found)}")
    if carbon_count > 0:
        reason_parts.append("simple inorganic/organic salt")
    else:
        reason_parts.append("inorganic compound")
    
    return True, " - ".join(reason_parts)