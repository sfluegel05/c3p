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
        'Ca', 'Mg', 'K', 'Na', 'Fe', 'Zn', 'Cu', 'Mn',  # Major minerals
        'Ba', 'Sr', 'La', 'Cs', 'Al', 'Sb'  # Less common but relevant
    }
    
    # Allowed inorganic anion components
    allowed_anions = {
        'Cl', 'F', 'Br', 'I',  # Halides
        'P', 'S', 'Si',  # Other common inorganic elements
        'N', 'O', 'H'  # For nitrates, phosphates, hydroxides
    }

    # Get all elements present
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

    # Check for essential mineral cations
    essential_cations_found = elements.intersection(essential_cations)
    if not essential_cations_found:
        return False, "No essential mineral cations found"

    # Count atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    # Define common inorganic anion patterns
    common_anion_patterns = [
        ('[O-]P([O-])([O-])=O', 'phosphate'),  # Phosphate
        ('[O-]S([O-])(=O)=O', 'sulfate'),      # Sulfate
        ('[O-]C([O-])=O', 'carbonate'),        # Carbonate
        ('[O-][N+]([O-])=O', 'nitrate'),       # Nitrate
        ('F', 'fluoride'),                      # Fluoride
        ('Cl', 'chloride'),                     # Chloride
        ('[OH-]', 'hydroxide'),                 # Hydroxide
        ('[O-][Si]', 'silicate')               # Silicate
    ]

    # Check for presence of common inorganic anions
    found_anions = []
    for pattern, name in common_anion_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
            found_anions.append(name)

    # Handle carbon-containing compounds
    if carbon_count > 0:
        # Allow only if it's a carbonate
        if not ('carbonate' in found_anions and carbon_count == found_anions.count('carbonate')):
            return False, "Contains organic components - not a mineral nutrient"

    # Check remaining elements are acceptable
    non_essential_elements = elements - essential_cations
    if not non_essential_elements.issubset(allowed_anions):
        unexpected_elements = non_essential_elements - allowed_anions
        return False, f"Contains unexpected elements: {unexpected_elements}"

    # Verify ionic nature or valid metal compound
    if not (has_cation or has_anion):
        # Check for direct metal-containing compounds
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in essential_cations:
                if not any(neigh.GetSymbol() in ['O', 'F', 'Cl', 'Br', 'I'] 
                          for neigh in atom.GetNeighbors()):
                    return False, "Invalid metal bonding pattern"

    # Build classification reason
    if found_anions:
        anion_str = " and ".join(set(found_anions))
        reason = f"Contains {anion_str} with essential elements: {', '.join(essential_cations_found)}"
    else:
        reason = f"Simple inorganic compound with essential elements: {', '.join(essential_cations_found)}"

    return True, reason