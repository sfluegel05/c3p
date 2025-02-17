"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral (geologically formed inorganic species)
Heuristics revised:
  1. Many minerals are salts made of disconnected (ionic) fragments.
  2. In a multi‐fragment species, fragments that are monoatomic (or water) are taken as inorganic.
  3. For multi‐atom fragments, if a fragment does not contain any carbon it is assumed inorganic.
  4. If a fragment does contain carbon, it is presumed organic unless it consists exclusively of C, O and H
     and has at least 6 carbon atoms (a “fatty acid” or simple carboxylate); or it is a small simple carboxylate.
  5. For single–fragment species we now require that:
         • If no carbon is present, then only small fragments (≤6 heavy atoms) are considered mineral.
         • If carbon is present, the species is taken as organic.
  6. In addition, when metals are present we “whitelist” only a set of common, geologically‐relevant metals.
     
Note: This is a very heuristic approach.
"""
from rdkit import Chem

# A small helper to test for water (H2O) fragments.
def is_water(frag):
    if frag.GetNumAtoms() == 3:
        nums = sorted([atom.GetAtomicNum() for atom in frag.GetAtoms()])
        if nums == [1, 1, 8]:
            return True
    return False

# Allowed metal atomic numbers (common rock-formers).
ALLOWED_METALS = {3, 11, 12, 13, 19, 20, 26, 27, 28, 29, 30, 51, 55, 56, 57}  # e.g. Li, Na, Mg, Al, K, Ca, Fe, Co, Ni, Cu, Zn, Sb, Cs, Ba, La

# For matching a simple carboxylate group.
CARBOXYLATE_SMARTS = "[CX3](=O)[O-]"
CARBOX_PAT = Chem.MolFromSmarts(CARBOXYLATE_SMARTS)

def classify_fragment(frag):
    """
    Returns a string "inorganic" or "organic" for a multi–atom fragment.
    Heuristics:
      - Monoatomic fragments should not reach here (handled separately).
      - Water is inorganic.
      - Fragments without any carbon are treated as inorganic if they do not contain any metals
        outside the allowed set.
      - For fragments that contain carbon:
            * If the fragment’s atoms (ignoring H) are exclusively from {C (6), O (8)} then:
                – if it has >= 6 C atoms, count it as inorganic (a fatty acid/carbonate–like ion).
                – if it is small (6 or fewer heavy atoms) and matches a carboxylate pattern, inorganic.
                – otherwise classify as organic.
            * If the fragment contains carbon plus any atom beyond {6,8,1} (and halogens may be allowed in organic):
                – If any metal is present that is not in ALLOWED_METALS, then flag as organic.
                – Otherwise, flag as organic.
      - Fragments that do not contain carbon are inorganic.
    """
    # Water fragments in multiatom form:
    if is_water(frag):
        return "inorganic"
    
    atoms = list(frag.GetAtoms())
    heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() != 1]
    heavy_count = len(heavy_atoms)
    elems = [atom.GetAtomicNum() for atom in heavy_atoms]

    # If no carbon is present, check metal whitelist.
    if 6 not in elems:
        # Even if metals are present, if any metal is not in our white list, mark as non-mineral.
        for atom in heavy_atoms:
            num = atom.GetAtomicNum()
            if num >= 3 and Chem.PeriodicTable.GetDefaultTable().IsMetal(num):
                if num not in ALLOWED_METALS:
                    return "organic"
        return "inorganic"
    
    # Fragment contains carbon. If all heavy atoms are only C and O, then it might be a simple anion.
    allowed_for_fatty = {6, 8}
    if set(elems).issubset(allowed_for_fatty):
        num_C = sum(1 for num in elems if num == 6)
        if num_C >= 6:
            return "inorganic"
        # Also, if it is small and matches a carboxylate pattern, treat as inorganic.
        if frag.HasSubstructMatch(CARBOX_PAT) and heavy_count <= 6:
            return "inorganic"
        return "organic"
    
    # Fragment contains carbon plus additional elements.
    # Now check if any metal present is not in the allowed list.
    for atom in heavy_atoms:
        num = atom.GetAtomicNum()
        if Chem.PeriodicTable.GetDefaultTable().IsMetal(num):
            if num not in ALLOWED_METALS:
                return "organic"
    # For safety, if carbon is present in a larger fragment, we consider it organic.
    return "organic"

def is_mineral(smiles: str):
    """
    Determines if a molecule should be considered a mineral based on its SMILES string
    using a heuristic that discriminates on the basis of disconnected fragments, 
    carbon content, and the type of metals present.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if molecule meets the minerals criteria, False otherwise.
       str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get disconnected fragments.
    frags = Chem.GetMolFrags(mol, asMols=True)
    num_frags = len(frags)
    
    # ----- Single Fragment Species -----
    if num_frags == 1:
        frag = frags[0]
        # Monoatomic species are not minerals.
        if frag.GetNumAtoms() < 2:
            return False, "Single atom species is not a mineral"
        # Count heavy atoms (non‐hydrogen)
        heavy_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1]
        heavy_count = len(heavy_atoms)
        # Check for carbon.
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms())
        if not has_carbon:
            # For single–fragment species lacking carbon, require a small heavy-atom count (e.g. ≤6)
            if heavy_count <= 6:
                return True, "Single fragment without carbon and small heavy-atom count indicates typical ionic mineral nature"
            else:
                return False, "Single fragment without carbon but heavy-atom count is high, suggesting a complex coordination not typical for minerals"
        else:
            return False, "Single fragment contains carbon, suggesting an organic species"
    
    # ----- Multiple Fragment Species -----
    # We separate fragments into monoatomic (or water) and multi–atom fragments.
    inorganic_frag_count = 0
    organic_frag_count = 0
    total_frag = num_frags
    for frag in frags:
        if frag.GetNumAtoms() == 1:
            inorganic_frag_count += 1
            continue
        if is_water(frag):
            inorganic_frag_count += 1
            continue
        # For multi–atom fragments, classify using our heuristic.
        flag = classify_fragment(frag)
        if flag == "inorganic":
            inorganic_frag_count += 1
        else:
            organic_frag_count += 1

    # Decide overall:
    # If at least half the fragments are inorganic we call it mineral.
    if inorganic_frag_count >= organic_frag_count:
        return True, f"Multiple fragments detected (count: {total_frag}) with predominantly inorganic character (inorganic: {inorganic_frag_count}, organic: {organic_frag_count})."
    else:
        return False, f"Multiple fragments detected (count: {total_frag}) with predominantly organic character (inorganic: {inorganic_frag_count}, organic: {organic_frag_count})."

# --- Example calls for testing (uncomment to run) ---
if __name__ == '__main__':
    # A set of test cases (as provided) with their SMILES strings.
    test_smiles = {
        "greigite": "[S--].[S--].[S--].[S--].[Fe++].[Fe+3].[Fe+3]",
        "calcium difluoride": "[F-].[F-].[Ca++]",
        "magnesium distearate": "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",
        "disodium tetraborate": "[Na+].[Na+].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2",
        "chalcopyrite": "[S--].[S--].[Fe++].[Cu++]",
        "caesium formate": "[Cs+].[H]C([O-])=O",
        "hexachloroiridate(2-)": "[Cl-][Ir--](Cl)(Cl)(Cl)(Cl)Cl",
        "cobalt-precorrin-6B": "[H][C@]12[C@H](CC(O)=O)[C@@](C)(CCC(O)=O)C3=[N+]1[Co--]14N5C(=CC6=[N+]1C(C[C@@]1(C)C(CC(O)=O)=C(CCC(O)=O)C(C3)=[N+]41)=C(CCC(O)=O)[C@]6(C)CC(O)=O)[C@@H](CCC(O)=O)[C@](C)(CC(O)=O)[C@]25C",
        "sodium tungstate": "[Na+].[Na+].[O-][W]([O-])(=O)=O",
        "antimony pentafluoride": "[Sb](F)(F)(F)(F)F",
        "monosodium urate": "[H]C1=NC2=NC(=O)NC(=O)N2C(=O)[O-].[Na+]",
        "tin(II) chloride (anhydrous)": "[Cl-].[Cl-].[Sn++]"
    }
    
    for name, smi in test_smiles.items():
        flag, reason = is_mineral(smi)
        print(f"{name}: {flag} – {reason}")