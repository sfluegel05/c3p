"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: CHEBI:27025 mineral nutrient
Definition: A mineral that is an inorganic nutrient which must be ingested and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or structural functions in the human body.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic compound containing essential mineral elements and mineral nutrient anions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """

    # Parse SMILES and split into ions
    mol_list = [Chem.MolFromSmiles(s) for s in smiles.split(".")]
    if any(mol is None for mol in mol_list):
        return False, "Invalid SMILES string"

    # Check for organic functional groups or long carbon chains
    organic_smarts = "[C:1]=[C:2]|[C:1]#[C:2]|[C:1]~[C:2]~[C:3]~[C:4]~[C:5]~[C:6]"
    organic_pattern = Chem.MolFromSmarts(organic_smarts)
    if any(mol.HasSubstructMatch(organic_pattern) for mol in mol_list):
        return False, "Contains organic functional groups or long carbon chains"

    # Check for essential mineral elements
    mineral_elements = [3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103]
    elements_present = set.union(*[set(atom.GetAtomicNum() for atom in mol.GetAtoms()) for mol in mol_list])
    mineral_elements_present = elements_present.intersection(mineral_elements)
    if not mineral_elements_present:
        return False, "Does not contain any essential mineral elements"

    # Check for mineral nutrient cations
    cation_patterns = [
        Chem.MolFromSmarts("[K+]"),   # potassium
        Chem.MolFromSmarts("[Ca++]"), # calcium
        Chem.MolFromSmarts("[Mg++]"), # magnesium
        Chem.MolFromSmarts("[Fe+2]"), # iron(2+)
        Chem.MolFromSmarts("[Fe+3]"), # iron(3+)
        Chem.MolFromSmarts("[Zn++]"), # zinc
        Chem.MolFromSmarts("[Cu+]"),  # copper(1+)
        Chem.MolFromSmarts("[Cu++]"), # copper(2+)
        Chem.MolFromSmarts("[Mn+2]"), # manganese(2+)
        Chem.MolFromSmarts("[Na+]"),  # sodium
        Chem.MolFromSmarts("[Cs+]"),  # cesium
        Chem.MolFromSmarts("[Ba++]"), # barium
        Chem.MolFromSmarts("[Al+3]"), # aluminium
        # Add more cation patterns as needed
    ]
    cation_matches = [any(mol.HasSubstructMatch(pattern) for mol in mol_list) for pattern in cation_patterns]
    if not any(cation_matches):
        return False, "Does not contain any common mineral nutrient cations"

    # Check for mineral nutrient anions
    anion_patterns = [
        Chem.MolFromSmarts("[O-]P([O-])([O-])=O"),  # phosphate
        Chem.MolFromSmarts("[O-]S([O-])(=O)=O"),    # sulfate
        Chem.MolFromSmarts("[Cl-]"),                # chloride
        Chem.MolFromSmarts("[F-]"),                 # fluoride
        Chem.MolFromSmarts("[O-]C([O-])=O"),        # carbonate
        Chem.MolFromSmarts("[O-][N+]([O-])=O"),     # nitrate
        Chem.MolFromSmarts("[O-][Si]([O-])([O-])[O-]"),  # silicate
        Chem.MolFromSmarts("[O-][B]([O-])([O-])"),      # borate
        Chem.MolFromSmarts("[O-H]"),                    # hydroxide
        Chem.MolFromSmarts("[O-2]"),                    # oxide
        # Add more anion patterns as needed
    ]
    anion_matches = [any(mol.HasSubstructMatch(pattern) for mol in mol_list) for pattern in anion_patterns]
    if not any(anion_matches):
        return False, "Does not contain any common mineral nutrient anions"

    # Check molecular weight and formal charge
    mol_wt_sum = sum(rdMolDescriptors.CalcExactMolWt(mol) for mol in mol_list)
    if mol_wt_sum < 100:
        return False, "Molecular weight too low for mineral nutrient"

    formal_charge_sum = sum(rdMolDescriptors.CalcFormalCharge(mol) for mol in mol_list)
    if formal_charge_sum != 0:
        return False, "Overall formal charge is not zero"

    return True, "Inorganic compound containing essential mineral elements, mineral nutrient cations, and mineral nutrient anions"