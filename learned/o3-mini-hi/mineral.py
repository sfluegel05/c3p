"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral
Heuristics for mineral classification:
  1. Many minerals are salts made of disconnected (ionic) fragments.
  2. For multi‐fragment species, fragments that are monoatomic (or water) are taken as inorganic.
  3. For multi‐atom fragments, if a fragment does not contain any carbon it is assumed inorganic.
  4. If a fragment does contain carbon, it is presumed organic unless it consists exclusively of C, O and H
     and has at least 6 carbon atoms (a “fatty acid” or simple carboxylate); or it is a small, simple carboxylate.
  5. For single–fragment species:
         • If no carbon is present, then only small fragments (≤6 heavy atoms) are considered mineral.
         • If carbon is present, the species is taken as organic.
  6. When metals are present, only a set of common, geologically‐relevant metals (white‐list) is allowed.
     
Note: This is a very heuristic approach.
"""

from rdkit import Chem

# Get periodic table (updated access method)
pt = Chem.GetPeriodicTable()

# A helper function to test for water (H2O) fragments.
def is_water(frag):
    # Water should have exactly three atoms: two H and one O.
    if frag.GetNumAtoms() == 3:
        atomic_nums = sorted([atom.GetAtomicNum() for atom in frag.GetAtoms()])
        if atomic_nums == [1, 1, 8]:
            return True
    return False

# Allowed metal atomic numbers (common rock-formers).
ALLOWED_METALS = {3, 11, 12, 13, 19, 20, 26, 27, 28, 29, 30, 51, 55, 56, 57}  # e.g. Li, Na, Mg, Al, K, Ca, Fe, Co, Ni, Cu, Zn, Sb, Cs, Ba, La

# For matching a simple carboxylate group.
CARBOXYLATE_SMARTS = "[CX3](=O)[O-]"
CARBOX_PAT = Chem.MolFromSmarts(CARBOXYLATE_SMARTS)

def classify_fragment(frag):
    """
    Classify a multi-atom fragment as 'inorganic' or 'organic' based on its atoms.
    Heuristics:
      - Water is inorganic.
      - Fragments without carbon are treated as inorganic if metals (if any) belong to the allowed list.
      - Fragments that contain carbon:
            * If all heavy atoms are only C and O: if there are at least 6 C, consider inorganic,
              or if small (≤6 heavy atoms) and matching a carboxylate, treat as inorganic.
            * Otherwise, classify as organic.
    """
    if is_water(frag):
        return "inorganic"
    
    heavy_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1]
    elems = [atom.GetAtomicNum() for atom in heavy_atoms]
    
    # If no carbon is present, check metals.
    if 6 not in elems:
        for atom in heavy_atoms:
            num = atom.GetAtomicNum()
            # If atom is a metal and not allowed, flag as organic.
            if num >= 3 and pt.IsMetal(num):
                if num not in ALLOWED_METALS:
                    return "organic"
        return "inorganic"
    
    # Fragment contains carbon.
    allowed_for_fatty = {6, 8}
    if set(elems).issubset(allowed_for_fatty):
        # Count carbons.
        num_C = sum(1 for num in elems if num == 6)
        if num_C >= 6:
            return "inorganic"  # Possibly a fatty acid or simple inorganic carboxylate
        # Also, if it is small and matches a carboxylate pattern, treat as inorganic.
        if frag.HasSubstructMatch(CARBOX_PAT) and len(heavy_atoms) <= 6:
            return "inorganic"
        return "organic"
    
    # Fragments with carbon plus other elements. Check if any metal is outside the allowed metal list.
    for atom in heavy_atoms:
        num = atom.GetAtomicNum()
        if pt.IsMetal(num):
            if num not in ALLOWED_METALS:
                return "organic"
    # Generally, if carbon is present in a more complex context, classify as organic.
    return "organic"

def is_mineral(smiles: str):
    """
    Detects if a molecule is a mineral based on a heuristic analysis of its fragments.
    
    Args:
       smiles (str): SMILES string for the molecule.
       
    Returns:
       bool: True if the molecule meets the mineral criteria, False otherwise.
       str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split molecule into disconnected fragments.
    frags = Chem.GetMolFrags(mol, asMols=True)
    num_frags = len(frags)
    
    # Single Fragment scenario.
    if num_frags == 1:
        frag = frags[0]
        if frag.GetNumAtoms() < 2:
            return False, "Single atom species are not considered minerals"
        heavy_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1]
        heavy_count = len(heavy_atoms)
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms())
        if not has_carbon:
            if heavy_count <= 6:
                return True, "Single fragment without carbon and small heavy-atom count suggests an ionic mineral."
            else:
                return False, "Single fragment without carbon but heavy-atom count is too high for simple minerals."
        else:
            return False, "Single fragment contains carbon, suggesting an organic species."
    
    # Multiple fragments scenario.
    inorganic_fragments = 0
    organic_fragments = 0
    for frag in frags:
        # Monoatomic fragments are considered inorganic.
        if frag.GetNumAtoms() == 1:
            inorganic_fragments += 1
            continue
        # Water is inorganic.
        if is_water(frag):
            inorganic_fragments += 1
            continue
        # Classify the multi–atom fragment with our heuristic.
        frag_type = classify_fragment(frag)
        if frag_type == "inorganic":
            inorganic_fragments += 1
        else:
            organic_fragments += 1
    
    total = len(frags)
    # Decide overall: if at least half the fragments are inorganic, classify as mineral.
    if inorganic_fragments >= organic_fragments:
        return True, f"Multiple fragments (total: {total}) with predominantly inorganic character (inorganic: {inorganic_fragments}, organic: {organic_fragments})."
    else:
        return False, f"Multiple fragments (total: {total}) with predominantly organic character (inorganic: {inorganic_fragments}, organic: {organic_fragments})."

# Optional: test cases to check the classification.
if __name__ == "__main__":
    test_smiles = {
        "greigite": "[S--].[S--].[S--].[S--].[Fe++].[Fe+3].[Fe+3]",
        "calcium difluoride": "[F-].[F-].[Ca++]",
        "magnesium distearate": "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",
        "disodium tetraborate": "[Na+].[Na+].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2",
        "chalcopyrite": "[S--].[S--].[Fe++].[Cu++]",
        "caesium formate": "[Cs+].[H]C([O-])=O",
        "antimony pentafluoride": "[Sb](F)(F)(F)(F)F",
        "calcium dihydroxide": "[OH-].[OH-].[Ca++]",
        "calcium carbonate": "[Ca+2].C(=O)([O-])[O-]",
        "disodium tetraborate octahydrate": "[Na+].[Na+].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H]",
        "kaolinite": "[OH-].[OH-].[OH-].[OH-].O=[Si]([O-])O[Si](=O)[O-].[Al+3].[Al+3]",
        "barium carbonate": "[Ba++].[O-]C([O-])=O"
    }

    for name, smi in test_smiles.items():
        result, explanation = is_mineral(smi)
        print(f"{name}: {result} — {explanation}")