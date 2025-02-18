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
            if atom.GetFormalCharge() > 0 and atom.GetIsMetal():
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

    return True, "Contains metal cation and anion with balanced charges"