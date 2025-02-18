"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer is defined as a mixture composed of macromolecules of different kinds,
which may be differentiated by composition, length, degree of branching etc.
This improved heuristic attempts to disregard small counterions and salts by requiring that
each fragment be “macromolecular” (e.g., above a certain molecular weight and heavy-atom count)
and by discarding fragments that are very small relative to the largest fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer – that is, a mixture of distinct
    macromolecular components.

    Heuristic improvements:
      1. Parse the SMILES and split it into disconnected fragments.
      2. Discard fragments that are too small to be considered macromolecules (require MW >= 200 Da
         and at least 10 heavy atoms).
      3. Discard fragments that, although meeting the basic criteria, are much smaller than the largest fragment.
         (For example, if one fragment is small enough to be a counterion, it will be ignored.)
      4. Among the remaining fragments, generate canonical SMILES and check whether there are at least
         two distinct entities.
         (If the mixture only contains multiple copies of the same fragment, it is more likely a salt than a polymer.)

    Args:
        smiles (str): SMILES string representing the substance.

    Returns:
        bool: True if classified as a polymer, False otherwise.
        str: Explanation for the classification decision.
    """

    # Parse the input; if invalid, return immediately.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Split the molecule into disconnected fragments.
    fragments = Chem.GetMolFrags(mol, asMols=True)

    # Allowed atomic numbers for common organic macromolecules
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}

    # Define a function that checks the “macromolecular” quality of a fragment.
    def is_macromolecular(frag):
        # Calculate molecular weight and count heavy atoms (exclude hydrogens)
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        heavy_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        # Discard if below thresholds
        if mw < 200 or heavy_atoms < 10:
            return False
        # Disqualify a fragment if it contains atoms outside our allowed set.
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() not in allowed_atomic_nums:
                return False
        return True

    # First, filter fragments using the basic macromolecular criteria.
    candidate_frags = [frag for frag in fragments if is_macromolecular(frag)]
    if not candidate_frags:
        return False, "No fragment qualifies as a macromolecule based on MW and heavy atom count"

    # Now, determine the maximum MW among these candidates.
    frag_mws = [rdMolDescriptors.CalcExactMolWt(frag) for frag in candidate_frags]
    max_mw = max(frag_mws)

    # Discard fragments that are very small relative to the largest fragment
    # Here, we demand that the fragment's MW be at least 60% of the max; otherwise it is likely a counterion.
    filtered_frags = []
    for frag, mw in zip(candidate_frags, frag_mws):
        if mw >= 0.6 * max_mw:
            filtered_frags.append(frag)

    # If fewer than 2 remain, it is not a mixture of distinct macromolecules.
    if len(filtered_frags) < 2:
        return False, ("Does not meet the criteria for a polymer "
                       "(fewer than two sufficiently large macromolecular fragments)")

    # Create a set of canonical SMILES from the filtered fragments to check for chemical distinctness.
    unique_frag_smiles = set(Chem.MolToSmiles(frag, canonical=True) for frag in filtered_frags)
    if len(unique_frag_smiles) >= 2:
        return True, ("Detected a mixture of at least two distinct macromolecular components "
                      "(each with MW >=200 Da, at least 10 heavy atoms, and not much smaller than the largest fragment)")
    else:
        # If the fragments are chemically identical, this is more likely a salt than a polymer.
        return False, ("Multiple large fragments were found but they appear chemically identical; "
                       "this is more likely a salt rather than a polymer mixture")

# Example usage (uncomment to test):
# test_smiles = "C=1(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N"
# result, reason = is_polymer(test_smiles)
# print(result, reason)