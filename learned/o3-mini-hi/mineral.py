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

# Define a helper function to check if a fragment is water (H2O).
def is_water(frag):
    # Water should have exactly three atoms: 2 H and 1 O.
    if frag.GetNumAtoms() == 3:
        atomic_nums = sorted([atom.GetAtomicNum() for atom in frag.GetAtoms()])
        if atomic_nums == [1, 1, 8]:
            return True
    return False

# Allowed metal atomic numbers (common rock-formers).
ALLOWED_METALS = {3, 11, 12, 13, 19, 20, 26, 27, 28, 29, 30, 51, 55, 56, 57}  # e.g. Li, Na, Mg, Al, K, Ca, Fe, Co, Ni, Cu, Zn, Sb, Cs, Ba, La

# Define a set of all metals (a simplified set covering alkali, alkaline earth, most transition metals, and some post-transition metals)
ALL_METALS = {
    # Alkali metals
    3, 11, 19, 37, 55, 87,
    # Alkaline earth metals
    4, 12, 20, 38, 56, 88,
    # Transition metals (first and second row, plus a few from later rows)
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    57, 72, 73, 74, 75, 76, 77, 78, 79, 80,
    # Some post-transition metals (if needed)
    13, 31, 49, 50, 81, 82, 83
}

# SMARTS pattern for a carboxylate group
CARBOXYLATE_SMARTS = "[CX3](=O)[O-]"
CARBOX_PAT = Chem.MolFromSmarts(CARBOXYLATE_SMARTS)

def classify_fragment(frag):
    """
    Classify a multi-atom fragment as 'inorganic' or 'organic' based on its atoms.
    """
    if is_water(frag):
        return "inorganic"
    
    heavy_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1]
    elems = [atom.GetAtomicNum() for atom in heavy_atoms]
    
    # If no carbon is present, check for metals.
    if 6 not in elems:
        for atom in heavy_atoms:
            num = atom.GetAtomicNum()
            # If an atom is a metal (i.e. in ALL_METALS) but not in the allowed list, flag the fragment.
            if num in ALL_METALS and num not in ALLOWED_METALS:
                return "organic"
        return "inorganic"
    
    # Fragment contains carbon.
    allowed_for_fatty = {6, 8}  # carbon and oxygen only.
    if set(elems).issubset(allowed_for_fatty):
        # Count number of carbons.
        num_C = sum(1 for num in elems if num == 6)
        if num_C >= 6:
            return "inorganic"  # Possibly a fatty acid or simple inorganic carboxylate.
        # If the fragment is small and matches a carboxylate group, then treat as inorganic.
        if frag.HasSubstructMatch(CARBOX_PAT) and len(heavy_atoms) <= 6:
            return "inorganic"
        return "organic"
    
    # If the fragment includes carbon along with other atoms,
    # check if any metal is present that is not in the allowed list.
    for atom in heavy_atoms:
        num = atom.GetAtomicNum()
        if num in ALL_METALS and num not in ALLOWED_METALS:
            return "organic"
    # Generally, if carbon is present in a more complex context, classify the fragment as organic.
    return "organic"

def is_mineral(smiles: str):
    """
    Determines whether a molecule is a mineral based on its SMILES string using heuristic rules.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if classified as a mineral, False otherwise.
        str: Explanation for the classification.
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
            return False, "Single atom species are not considered minerals."
        heavy_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1]
        heavy_count = len(heavy_atoms)
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms())
        if not has_carbon:
            if heavy_count <= 6:
                return True, "Single fragment without carbon and small heavy-atom count suggests an ionic mineral."
            else:
                return False, "Single fragment without carbon but heavy-atom count is too high for a simple mineral."
        else:
            return False, "Single fragment contains carbon, suggesting an organic species."
    
    # Multiple fragment scenario.
    inorganic_fragments = 0
    organic_fragments = 0
    for frag in frags:
        # Monoatomic fragments are deemed inorganic.
        if frag.GetNumAtoms() == 1:
            inorganic_fragments += 1
            continue
        # Water fragments are inorganic.
        if is_water(frag):
            inorganic_fragments += 1
            continue
        # Classify the fragment based on its atom composition.
        frag_type = classify_fragment(frag)
        if frag_type == "inorganic":
            inorganic_fragments += 1
        else:
            organic_fragments += 1
    
    total = len(frags)
    # Decide overall: if at least as many fragments are inorganic compared to organic, 
    # we classify the molecule as a mineral.
    if inorganic_fragments >= organic_fragments:
        return True, f"Multiple fragments (total: {total}) with predominantly inorganic character (inorganic: {inorganic_fragments}, organic: {organic_fragments})."
    else:
        return False, f"Multiple fragments (total: {total}) with predominantly organic character (inorganic: {inorganic_fragments}, organic: {organic_fragments})."

# Optional test cases.
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