"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import Mol

# Expanded list of metal atomic numbers (including those from examples)
METAL_ATOMIC_NUMBERS = {
    # Alkali metals
    3, 11, 19, 37, 55,
    # Alkaline earth
    4, 12, 20, 38, 56,
    # Transition metals
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    72, 73, 74, 75, 76, 77, 78, 79, 80,
    # Post-transition metals
    13, 31, 49, 50, 51, 81, 82, 83,
    # Lanthanides
    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71
}

# Predefined anion patterns
CARBONYLATE = Chem.MolFromSmarts("[CX3](=O)[O-]")
CARBONATE = Chem.MolFromSmarts("[O-]C(=O)[O-]")

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic salt containing at least one metal cation and one anion,
    or a single-component ionic compound with a metal bonded to electronegative atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    
    # Case 1: Single-component ionic compound (e.g., SbF5, LaCl3)
    if len(frags) == 1:
        # Check for metal with positive charge
        metals = [atom for atom in mol.GetAtoms() 
                  if atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS 
                  and atom.GetFormalCharge() > 0]
        if not metals:
            return False, "No metal cation in single-component structure"
        
        # Check all metal bonds are to F, Cl, O, S
        for metal in metals:
            for nbr in metal.GetNeighbors():
                if nbr.GetAtomicNum() not in {9, 17, 8, 16}:  # F, Cl, O, S
                    return False, "Metal bonded to non-electronegative atoms"
        
        # Check overall charge balance
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        if total_charge != 0:
            return False, "Charge imbalance in single-component structure"
        
        # Check for any carbon atoms
        if any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
            return False, "Single-component structure contains carbon"
        
        return True, "Single-component ionic compound with metal and electronegative anions"
    
    # Case 2: Multi-component salt
    total_charge = 0
    has_metal_cation = False
    has_anion = False
    
    for frag in frags:
        frag_charge = Chem.GetFormalCharge(frag)
        total_charge += frag_charge
        
        # Check for metal cations (can be multi-atom like NH4+ but we focus on single-atom)
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS and atom.GetFormalCharge() > 0:
                has_metal_cation = True
        
        # Check anions (negative charge fragments)
        if frag_charge < 0:
            has_anion = True
            # Check if anion contains carbon that's not in carbonate/carboxylate
            if any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms()):
                # Check for carboxylate or carbonate patterns
                has_carboxylate = frag.HasSubstructMatch(CARBONYLATE)
                has_carbonate = frag.HasSubstructMatch(CARBONATE)
                if not (has_carboxylate or has_carbonate):
                    return False, "Anion contains non-allowed organic components"
    
    # Check charge balance
    if total_charge != 0:
        return False, f"Charge imbalance ({total_charge})"
    
    if not has_metal_cation:
        return False, "No metal cation detected"
    
    if not has_anion:
        return False, "No anion detected"
    
    # Check for organic components outside of allowed anions
    for frag in frags:
        if Chem.GetFormalCharge(frag) >= 0:  # Check cations and neutral fragments
            if any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms()):
                return False, "Organic components in cation/neutral fragments"
    
    return True, "Contains metal cation and anion with balanced charges"