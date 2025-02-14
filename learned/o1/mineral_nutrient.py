"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.

    A mineral nutrient is a mineral that is an inorganic nutrient which must be ingested 
    and absorbed in adequate amounts to satisfy a wide range of essential metabolic
    and/or structural functions in the human body.

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

    # List of essential mineral elements
    essential_minerals = {
        "Ca",  # Calcium
        "P",   # Phosphorus
        "K",   # Potassium
        "Na",  # Sodium
        "Cl",  # Chlorine
        "Mg",  # Magnesium
        "S",   # Sulfur
        "Fe",  # Iron
        "Mn",  # Manganese
        "Cu",  # Copper
        "I",   # Iodine
        "Zn",  # Zinc
        "Se",  # Selenium
        "Mo",  # Molybdenum
        "Cr",  # Chromium
        "Co",  # Cobalt
        "F",   # Fluorine
        "Ba",  # Barium
        "Sb",  # Antimony
        "La",  # Lanthanum
        "Cs",  # Cesium
        "Pd",  # Palladium
        "Al",  # Aluminium
        "Si"   # Silicon
    }

    # Split molecule into disconnected fragments (ions)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    has_essential_mineral_cation = False

    for frag in frags:
        atoms = frag.GetAtoms()
        symbols = {atom.GetSymbol() for atom in atoms}
        atomic_nums = {atom.GetAtomicNum() for atom in atoms}
        has_carbon = 6 in atomic_nums

        # Check if fragment contains essential mineral elements
        essential_mineral_in_frag = symbols.intersection(essential_minerals)
        if essential_mineral_in_frag:
            if has_carbon:
                # Essential mineral is part of an organic fragment
                continue
            else:
                # Fragment is a cation containing an essential mineral not bonded to carbon
                has_essential_mineral_cation = True
                # No need to check further fragments
                break

    if not has_essential_mineral_cation:
        return False, "Does not contain essential mineral cation not bonded to carbon"

    return True, "Contains essential mineral cation not bonded to carbon"