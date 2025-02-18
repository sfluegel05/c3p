"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    
    A mineral is typically an inorganic compound that is crystalline in nature and formed through geological processes.
    For this function, minerals are characterized by:
    - Absence of organic functional groups (no C-H or C-C bonds).
    - Presence of metal atoms.
    - May contain common inorganic anions (e.g., sulfate, phosphate).
    - May include water molecules (hydrates).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a mineral, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of metal atomic numbers (simplified for common metals in minerals)
    metal_atomic_numbers = [
        3,  # Lithium
        4,  # Beryllium
        11, # Sodium
        12, # Magnesium
        13, # Aluminium
        19, # Potassium
        20, # Calcium
        21, # Scandium
        22, # Titanium
        23, # Vanadium
        24, # Chromium
        25, # Manganese
        26, # Iron
        27, # Cobalt
        28, # Nickel
        29, # Copper
        30, # Zinc
        31, # Gallium
        37, # Rubidium
        38, # Strontium
        39, # Yttrium
        40, # Zirconium
        47, # Silver
        48, # Cadmium
        49, # Indium
        55, # Caesium
        56, # Barium
        57, # Lanthanum
        72, # Hafnium
        79, # Gold
        80, # Mercury
        82, # Lead
        83, # Bismuth
    ]
    
    # Check for presence of metals
    metals_in_molecule = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metal_atomic_numbers:
            metals_in_molecule.add(atom.GetSymbol())
    if not metals_in_molecule:
        return False, "No metal atoms found in the molecule"
    
    # Check for presence of carbon
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    
    # Check for C-H and C-C bonds (indicative of organic molecules)
    has_c_h_bonds = False
    has_c_c_bonds = False
    if has_carbon:
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Check for C-H bond
            if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 1) or \
               (a1.GetAtomicNum() == 1 and a2.GetAtomicNum() == 6):
                has_c_h_bonds = True
                break
            # Check for C-C bond
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                has_c_c_bonds = True
                break
    if has_c_h_bonds or has_c_c_bonds:
        return False, "Contains C-H or C-C bonds, indicative of organic molecule"
    
    # Check for common inorganic anions (e.g., sulfate, phosphate)
    inorganic_anions = {
        "sulfate": Chem.MolFromSmarts("S(=O)(=O)([O-])[O-]"),
        "phosphate": Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"),
        "carbonate": Chem.MolFromSmarts("C(=O)([O-])[O-]"),
        "nitrate": Chem.MolFromSmarts("N(=O)([O-])[O-]"),
        "chloride": Chem.MolFromSmarts("[Cl-]"),
        "fluoride": Chem.MolFromSmarts("[F-]"),
        "hydroxide": Chem.MolFromSmarts("[OH-]"),
        "oxide": Chem.MolFromSmarts("[O-2]")
    }
    anions_found = []
    for name, pattern in inorganic_anions.items():
        if mol.HasSubstructMatch(pattern):
            anions_found.append(name)
    
    # Check for hydrates (water molecules)
    water_pattern = Chem.MolFromSmarts("O")
    has_water = mol.HasSubstructMatch(water_pattern)
    
    # If all checks pass, classify as mineral
    reason = f"Classified as mineral: contains metals ({', '.join(metals_in_molecule)})"
    if anions_found:
        reason += f" and common inorganic anions ({', '.join(anions_found)})"
    if has_water:
        reason += "; includes water molecules (hydrate)"
    return True, reason