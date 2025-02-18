"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem

# Common metal elements in mineral nutrients (atomic numbers)
METAL_ATOMIC_NUMBERS = {
    3, 11, 19, 37, 55,  # Alkali metals: Li, Na, K, Rb, Cs
    4, 12, 20, 38, 56,  # Alkaline earth: Be, Mg, Ca, Sr, Ba
    13, 31, 49, 81,     # Group 13: Al, Ga, In, Tl
    26, 24, 25, 27, 28, 29, 30,  # Transition metals: Fe, Cr, Mn, Co, Ni, Cu, Zn
    57, 58, 59, 60      # Lanthanides: La, Ce, Pr, Nd
}

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic salt containing at least one metal cation and one anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Split into disconnected components (ions)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) < 2:
        return False, "Not a salt (single component)"

    total_charge = 0
    has_metal_cation = False
    has_anion = False

    for frag in frags:
        charge = Chem.GetFormalCharge(frag)
        total_charge += charge

        # Check for metal cations (single atom with positive charge)
        if frag.GetNumAtoms() == 1:
            atom = frag.GetAtomWithIdx(0)
            if atom.GetFormalCharge() > 0 and atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS:
                has_metal_cation = True
        
        # Check for anions (any component with negative charge)
        if charge < 0:
            has_anion = True

    # Check charge balance and required components
    if total_charge != 0:
        return False, f"Charge imbalance ({total_charge})"

    if not has_metal_cation:
        return False, "No metal cation detected"

    if not has_anion:
        return False, "No anion detected"

    # Basic check for organic components (crude but helps exclude organics)
    if any(atom.GetAtomicNum() == 6 for frag in frags for atom in frag.GetAtoms()):
        # Allow carbonate/bicarbonate and simple carboxylates (like acetate)
        # Check if any carbon has non-oxygen neighbors (indicates organic structure)
        for frag in frags:
            for atom in frag.GetAtoms():
                if atom.GetAtomicNum() == 6:
                    if any(nbr.GetAtomicNum() not in {8, 16, 15} for nbr in atom.GetNeighbors()):
                        return False, "Contains organic components"

    return True, "Contains metal cation and anion with balanced charges"