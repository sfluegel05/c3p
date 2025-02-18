"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer is defined as a mixture composed of macromolecules of different kinds,
which may be differentiated by composition, length, degree of branching etc.
This program implements an improved heuristic for classifying a SMILES string as a polymer.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer—that is, a mixture of distinct macromolecules.
    
    Heuristic improvements:
      1. Parse the SMILES and split it into disconnected fragments.
      2. Discard fragments that are too small (MW < 150 Da or fewer than 5 heavy atoms).
      3. Discard fragments that contain metals or elements not commonly found in organic macromolecules.
         We define an “allowed element” set as common organic elements: H, C, N, O, F, P, S, Cl, Br, I.
      4. From the remaining fragments, check if there are at least two that are chemically distinct.
         (If the mixture only contains multiple copies of the same large fragment, it is more likely a salt.)
    
    Args:
        smiles (str): SMILES string representing the substance.
        
    Returns:
        bool: True if classified as a polymer, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the given SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split the molecule into disconnected fragments.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    # Define allowed atomic numbers for organic macromolecules:
    # (includes H, C, N, O, F, P, S, Cl, Br, I)
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    
    # Function to check if a fragment is large enough and organic.
    def is_large_fragment(frag):
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        # Count heavy atoms (exclude hydrogens)
        heavy_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        if mw < 150 or heavy_atoms < 5:
            return False
        # Check that every atom is one of the allowed (organic) elements.
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() not in allowed_atomic_nums:
                # If a fragment contains a metal or uncommon element, skip it.
                return False
        return True

    # Filter fragments to keep only those that meet our threshold.
    large_frags = [frag for frag in fragments if is_large_fragment(frag)]
    
    if len(large_frags) < 2:
        return False, ("Does not meet the criteria for a polymer "
                       "(fewer than two distinct macromolecular fragments)")
    
    # Generate canonical SMILES for each large fragment.
    frag_smiles = [Chem.MolToSmiles(frag, canonical=True) for frag in large_frags]
    unique_fragments = set(frag_smiles)
    
    if len(unique_fragments) >= 2:
        return True, ("Detected a mixture of at least two distinct macromolecular components "
                      "(each with molecular weight >=150 Da, at least 5 heavy atoms, "
                      "and composed of common organic elements)")
    else:
        # Optionally: even if only one unique fragment is found, check if the copies vary in size.
        # For now, we consider multiple identical fragments as a salt (or non-polymer) mixture.
        return False, ("Multiple large fragments were found but they appear chemically identical; "
                       "this is more likely a salt rather than a polymer mixture")

# Example usage (uncomment the following lines to test):
# test_smiles = "C=1(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N"
# result, reason = is_polymer(test_smiles)
# print(result, reason)