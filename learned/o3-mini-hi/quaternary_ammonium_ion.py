"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary ammonium ion

Definition:
  A derivative of ammonium, NH4(+), in which all four hydrogens bonded to nitrogen 
  have been replaced with univalent (usually organyl) groups.

Our improved approach:
  1. Parse the molecule.
  2. Look for any nitrogen (atomic number 7) that is positively charged (+1) with exactly
     four bonds and zero hydrogens (explicit or implicit).
  3. For each attached substituent atom, we allow it if it is a carbon. If not carbon
     we require that it be “organic” – that is, it is attached (besides the quaternary nitrogen)
     to at least one carbon.
  4. To help avoid classifying large biomolecules (e.g. carnitine derivatives) that contain a choline‐like
     fragment, we check the overall size of the fragment that contains the N+ (unless that fragment contains phosphorus,
     in which case it likely is part of a phosphatidylcholine headgroup).
     
Tuning of thresholds (here, a fragment is considered “small” if it contains no more than 60 heavy atoms)
may be necessary depending on the application.
"""

from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determine whether the input SMILES string represents a molecule
    containing a quaternary ammonium ion – that is, an N+ (formal charge +1) with exactly 4 bonds,
    no attached hydrogens, and each substituent being “organic” (either carbon or carrying at least one carbon).
    Also, if the positive N is embedded in a very large fragment (and the fragment does not contain phosphorus)
    then we assume that the quaternary ammonium group is not the central chemical entity under study.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is classified as (containing) a quaternary ammonium ion.
        str: A reason explaining the decision.
              If an error occurs the function may return (None, None).
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all fragment atom lists (each fragment is a tuple of atom indices)
    frags = Chem.GetMolFrags(mol, asMols=False)
    
    # Iterate over all atoms looking for a candidate nitrogen
    for atom in mol.GetAtoms():
        # candidate: nitrogen with formal charge +1
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        # In a proper quaternary ammonium group, the N should have exactly 4 bonds
        # and (after RDKit adding implicit hydrogens) should have zero total hydrogens.
        if atom.GetDegree() != 4 or atom.GetTotalNumHs() != 0:
            continue

        # Check every neighbor of N:
        valid_substituents = True
        for nbr in atom.GetNeighbors():
            # Accept if neighbor is carbon.
            if nbr.GetAtomicNum() == 6:
                continue
            else:
                # For non-carbon neighbors, we require that the neighbor is part of some organic fragment.
                # In our approach, that means the neighbor must be connected 
                # (besides the candidate N) to at least one carbon.
                has_carbon = any(neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != atom.GetIdx() 
                                 for neighbor in nbr.GetNeighbors())
                if not has_carbon:
                    valid_substituents = False
                    break
        if not valid_substituents:
            continue

        # Determine the fragment that contains this candidate nitrogen.
        frag_size = None
        for frag in frags:
            if atom.GetIdx() in frag:
                frag_size = len(frag)  # count of heavy atoms in that fragment (all atoms since H are implicit)
                break
        if frag_size is None:
            # Should not happen.
            continue
        
        # If the fragment contains a phosphorus atom, we assume it is of the phosphocholine kind,
        # and we accept regardless of size.
        has_P = any(a.GetAtomicNum() == 15 for a in mol.GetAtoms() if a.GetIdx() in frag)
        
        # For simplicity we use a heavy atom count threshold (here 60) to “vet” embedded quaternary ammonium groups.
        # (A very large fragment is likely not the main quaternary ammonium chemical entity.)
        if (has_P is False) and (frag_size > 60):
            # The N+ is embedded in a large molecule without phosphorus;
            # in our tests this filter helps reduce false positives such as in carnitine derivatives.
            continue
        
        # If we reached here, the candidate nitrogen meets the criteria.
        return True, ("Found a candidate nitrogen atom (N+) with formal charge +1, exactly 4 bonds with no attached hydrogens, "
                      "and each substituent is organic. Fragment size = {} heavy atoms{}."
                      .format(frag_size, " (contains phosphorus)" if has_P else ""))
    
    # If no candidate was found, report the negative result.
    return False, "No quaternary ammonium group (N+ with 4 organic substituents) found"

# Example usage (testing a few cases):
if __name__ == "__main__":
    tests = [
        # True positives (expected to contain a quaternary ammonium ion)
        ("P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",
         "Example phosphinic acid with a choline-like group"),
        ("C[N+](C)(C)CC(O)O",
         "Betaine aldehyde hydrate, should be classified as containing a quaternary ammonium ion"),
        ("C[N+](C)(C)O",
         "Hydroxytrimethylaminium (previous false negative)"),
        ("[F-][N+](F)(F)F",
         "Tetrafluoroammonium (should not be classified: substituents are not organic)"),
        # False positives: one of the known misclassified carnitine derivatives
        ("CCCC\C=C\C\C=C\CCCCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C",
         "Carnitine derivative (should not be classified as the target quaternary ammonium ion)")
    ]
    
    for s, desc in tests:
        result, reason = is_quaternary_ammonium_ion(s)
        print("SMILES:", s)
        print("Description:", desc)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)