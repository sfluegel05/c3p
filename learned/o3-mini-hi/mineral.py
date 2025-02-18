"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral (inorganic, often ionic/geologically formed species)
Revised heuristic criteria:
  1. A mineral is typically an ionic salt (often hydrated) with more than one disconnected fragment.
  2. Even single‐fragment species may be mineral if they contain elements not common to organic chemistry.
  3. For multi–fragment species, if many (or the only multi–atom) fragments are “organic” (e.g. typical carboxylic acids or decorated organic groups)
     then the molecule is likely not a mineral. However, small simple inorganic anions (for example, carbonate 
     or simple long‐chain fatty acid anions having only C and O) are allowed.
     
Heuristics used below:
  – If the molecule consists only of disconnected single–atom “ions”, then (if two or more different elements are present)
    it is classified as mineral (calcium difluoride is an example), but a bare monoatomic ion is rejected.
  – In multi–fragment cases, for any fragment with more than one atom:
         • if the fragment contains any atom that is not in the “organic set” 
           (organic_set = {H, C, N, O, F, P, S, Cl, Br, I}) then it is presumed inorganic.
         • if the fragment does use only organic atoms but contains carbon then we look at its “size” and its
           “decorations”. In particular, if it is very small (or very decorated) then it is considered organic.
           Here we use two heuristics:
             (a) If a fragment contains carbon and its heavy–atom count is 6 or fewer AND it matches a simple 
                 carboxylate ([CX3](=O)[O-]) pattern then we let it count as inorganic.
             (b) Otherwise, if a fragment contains only C and O (and H) and has at least 6 carbons then we treat it as an
                 “inorganic‐type” fatty acid anion.
             (c) Otherwise it is taken to be organic.
  – For single–fragment species we require that (a) the molecule is not monoatomic and (b) that if it contains only atoms from
    the organic set then it is not considered mineral.
    
Note: This is only one (admittedly “messy” and heuristic) approach.
"""

from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is considered a mineral based on its SMILES string.
    
    The heuristic implemented here works roughly as follows:
      1. Parse the SMILES string.
      2. Get the disconnected fragments.
      3. For a molecule made of a single fragment:
            - Reject if it is a single (monoatomic) ion.
            - Otherwise, if it contains any atom outside a small “organic set” then consider it mineral.
      4. For a molecule made of multiple fragments:
            - If all fragments are monoatomic, then if there is more than one distinct element, classify as mineral.
            - For multi–atom fragments, check each fragment:
                  • if any atom is not in the organic set (organic_set = {H, C, N, O, F, P, S, Cl, Br, I}),
                    the fragment is flagged as inorganic.
                  • if the fragment uses only organic atoms but does contain C, then:
                        – if it matches a simple carboxylate pattern and is small (≤6 heavy atoms), 
                          count it as inorganic (e.g. formate, acetate, propionate).
                        – or if it consists solely of C, O (and H) and has at least 6 carbon atoms,
                          count it as inorganic (i.e. a simple fatty acid anion).
                        – otherwise, flag it as organic.
            - If more than half of the multi–atom fragments are flagged as organic, then the salt is likely not a mineral.
              Otherwise, it is.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule meets our criteria for a mineral, False otherwise.
       str: A reason describing the classification.
    """
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the set of elements common in organic molecules.
    organic_set = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}

    # For matching simple carboxylate groups
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pat = Chem.MolFromSmarts(carboxylate_smarts)

    # A helper to decide if a fragment is “water”
    def is_water(frag):
        # If the fragment has 3 atoms and exactly one oxygen and two hydrogens.
        if frag.GetNumAtoms() == 3:
            nums = sorted([atom.GetAtomicNum() for atom in frag.GetAtoms()])
            if nums == [1, 1, 8]:
                return True
        return False

    # A helper for multi-atom fragments containing C:
    def fragment_organic_flag(frag):
        """
        For fragments with more than one atom:
         - If any atom is not in the organic set, then treat as inorganic.
         - Else, if the fragment contains carbon then further decide:
             • If it is a simple carboxylate (matches pattern) and is small (≤6 heavy atoms) then treat inorganic.
             • Else if it consists solely of C, O, and H and has >= 6 carbon atoms then treat as inorganic.
             • Otherwise, flag as organic.
         - If no carbon is present, treat it as inorganic.
        Returns:
          True if the fragment is flagged as “organic” (i.e. not mineral-like), False otherwise.
        """
        # For water, force inorganic.
        if is_water(frag):
            return False

        atoms = list(frag.GetAtoms())
        # If any atom is not in the organic set, then fragment is inorganic.
        for atom in atoms:
            if atom.GetAtomicNum() not in organic_set:
                return False
        # Now, all atoms are drawn from the organic set.
        # If no carbon, then treat as inorganic.
        if not any(atom.GetAtomicNum() == 6 for atom in atoms):
            return False

        # Count heavy atoms (not hydrogen)
        heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() != 1]
        num_heavy = len(heavy_atoms)
        # Count number of carbon atoms
        num_C = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        # Count oxygens
        num_O = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)

        # If the fragment is a simple carboxylate (matches pattern) and is small (≤6 heavy atoms), treat as inorganic.
        if frag.HasSubstructMatch(carboxylate_pat) and num_heavy <= 6:
            return False

        # If the fragment is composed solely of C, O, and H and has at least 6 carbon atoms, treat as inorganic.
        atom_nums = set(atom.GetAtomicNum() for atom in atoms if atom.GetAtomicNum() != 1)
        if atom_nums.issubset({6, 8}) and num_C >= 6:
            return False

        # Otherwise, flag as organic.
        return True

    # Get the disconnected fragments (each as a separate molecule)
    frags = Chem.GetMolFrags(mol, asMols=True)
    frag_count = len(frags)

    # ----- Single fragment case -----
    if frag_count == 1:
        # If there is only one fragment and it is monoatomic, reject.
        if mol.GetNumAtoms() < 2:
            return False, "Single atom species is not a mineral"
        # If any atom is outside the organic set, classify as mineral.
        if any(atom.GetAtomicNum() not in organic_set for atom in mol.GetAtoms()):
            return True, "Single fragment with non‐organic atoms indicating possible mineral nature"
        else:
            return False, "Single fragment composed exclusively of typical organic elements"

    # ----- Multiple fragment case -----
    # If all fragments are single atoms, then require at least two different elements.
    single_atom_frags = [frag for frag in frags if frag.GetNumAtoms() == 1]
    multi_atom_frags = [frag for frag in frags if frag.GetNumAtoms() > 1]

    if len(frags) == len(single_atom_frags):
        distinct_elems = set()
        for frag in single_atom_frags:
            for atom in frag.GetAtoms():
                distinct_elems.add(atom.GetAtomicNum())
        if len(distinct_elems) > 1:
            return True, f"Multiple monoatomic fragments with {len(distinct_elems)} distinct elements typical for ionic mineral compounds."
        else:
            return False, "Multiple fragments but all are the same monoatomic ion"

    # For multi–fragment species that include at least one multi–atom fragment:
    organic_fragments = 0
    total_multi = 0
    for frag in multi_atom_frags:
        total_multi += 1
        # Determine if this multi–atom fragment is “organic” versus inorganic–like.
        if fragment_organic_flag(frag):
            organic_fragments += 1

    # Now decide:
    # If more than half of the multi-atom fragments are flagged as organic, then likely not a mineral.
    if total_multi > 0 and organic_fragments / total_multi > 0.5:
        return False, f"Majority ({organic_fragments} of {total_multi}) multi–atom fragments appear organic."
    else:
        return True, f"Multiple fragments detected (count: {frag_count}) with predominantly inorganic character typical for ionic mineral compounds."

# --- Example calls for testing – uncomment to run ---
if __name__ == '__main__':
    test_smiles = {
        "greigite": "[S--].[S--].[S--].[S--].[Fe++].[Fe+3].[Fe+3]",
        "calcium difluoride": "[F-].[F-].[Ca++]",
        "magnesium distearate": "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",
        "disodium tetraborate": "[Na+].[Na+].[H]OB1O[B-]2(O[H])OB(O[H])O[B-](O[H])(O1)O2",
        "chalcopyrite": "[S--].[S--].[Fe++].[Cu++]",
        "caesium formate": "[Cs+].[H]C([O-])=O",
        "potassium L-tartrate": "[K+].[K+].O[C@H]([C@@H](O)C([O-])=O)C([O-])=O"
    }
    for name, smi in test_smiles.items():
        is_min, reason = is_mineral(smi)
        print(f"{name}: {is_min} – {reason}")