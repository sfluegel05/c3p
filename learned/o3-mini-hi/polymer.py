"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer is defined as a mixture composed of macromolecules of different kinds,
which may be differentiated by composition, length, degree of branching etc.
This heuristic attempts further improvements by (a) requiring fragments to be big enough
(in terms of molecular weight, heavy-atom count, and a minimum fraction of the total mass)
and (b) comparing the chemical fingerprints of the candidate fragments to ensure they are sufficiently different.
If the improvements are insufficient, the function may return (None, None).
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors, AllChem

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer – that is, a mixture of distinct
    macromolecular components.
    
    Heuristic improvements:
      1. Parse the SMILES and split it into disconnected fragments.
      2. For each fragment compute its molecular weight (MW) and heavy atom count (HAC).
      3. Consider only fragments that meet the macromolecular criteria:
             - MW >= 300 Da
             - Heavy atom count >= 15
         and that contribute at least 15% of the total molecular weight.
      4. If there are at least two remaining fragments, compute Morgan fingerprints (radius=2, nBits=2048)
         and calculate pairwise Tanimoto similarity.
      5. If at least one pair of fragments has Tanimoto similarity below 0.85 (i.e. they are chemically distinct)
         then the SMILES is classified as polymer.
      6. Otherwise, if the fragments are all very similar (or if too few large fragments exist),
         the compound is not classified as a polymer.
    
    Args:
        smiles (str): SMILES string representing the substance.
    
    Returns:
        bool: True if classified as a polymer, False otherwise.
        str: Explanation for the classification decision.
    """
    # Try to parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split into disconnected fragments
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if not fragments:
        return False, "No fragments found in the molecule"
    
    # Compute the total molecular weight across ALL fragments.
    total_mw = sum(rdMolDescriptors.CalcExactMolWt(frag) for frag in fragments)
    
    # List to collect candidate fragments which are “macromolecular”
    candidate_frags = []
    candidate_data = []
    for frag in fragments:
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        # Count heavy atoms (non-hydrogen)
        heavy_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        # Consider fragment only if it meets our thresholds
        if mw >= 300 and heavy_atoms >= 15:
            # Also require that the fragment contributes at least 15% of the total MW
            if mw / total_mw >= 0.15:
                candidate_frags.append(frag)
                candidate_data.append({'mol': frag, 'mw': mw, 'heavy_atoms': heavy_atoms})
    
    if len(candidate_frags) < 2:
        return False, ("Does not meet the criteria for a polymer "
                       "(fewer than two sufficiently large macromolecular fragments)")
    
    # Now get canonical SMILES for distinctness (as a fallback)
    unique_smiles = set(Chem.MolToSmiles(frag, canonical=True) for frag in candidate_frags)
    if len(unique_smiles) < 2:
        return False, ("Multiple large fragments were found but they appear chemically identical; "
                       "this is more likely a salt rather than a polymer mixture")
    
    # Compute fingerprints for each candidate fragment
    fps = [AllChem.GetMorganFingerprintAsBitVect(frag, 2, nBits=2048) for frag in candidate_frags]
    
    # Check pairwise Tanimoto similarity; if at least one pair is sufficiently dissimilar, accept as polymer.
    n = len(fps)
    found_different = False
    for i in range(n):
        for j in range(i+1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            # If similarity is less than threshold (i.e., they share less than ~85% bits), then they are distinct.
            if sim < 0.85:
                found_different = True
                break
        if found_different:
            break

    if found_different:
        return True, ("Detected a mixture of at least two distinct macromolecular components "
                      "(each satisfying MW >=300 Da, heavy atoms >=15, and representing >=15% of total mass)")
    else:
        return False, ("All sufficiently large fragments appear very similar; likely a salt rather than a polymer mixture")
    
# Example usage (uncomment to test):
# test_smiles = "C=1C=C2C(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N"
# result, reason = is_polymer(test_smiles)
# print(result, reason)