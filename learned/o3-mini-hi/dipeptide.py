"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected via peptide linkages.

Improvements over previous attempt:
  - Instead of rejecting molecules with more than one candidate backbone bond we pick the candidate cleavage
    that yields the two fragments with the most balanced heavy atom counts.
  - We narrow the allowed heavy atom count for each fragment (set here between 4 and 35).
  - We continue to check for cyclic dipeptides (diketopiperazines) as an alternate classification.
  
The strategy:
  1. Parse the molecule and add hydrogens to make the connectivity more explicit.
  2. Check overall size; if too many heavy atoms, it is unlikely a dipeptide.
  3. Find all candidate amide bonds using a SMARTS pattern ([C](=O)[N]).
     For each match, if the bond is not aromatic, not in a ring, break the bond and check that exactly two fragments appear.
  4. For each candidate, count heavy atoms (atomic number > 1) in each fragment.
     Only consider candidates where both fragments have between 4 and 35 heavy atoms.
  5. If more than one candidate qualifies (which can happen when side chain amides also match), choose the candidate where the
     difference in heavy atom count between fragments is minimized.
  6. If one candidate is selected, return True with the explanation that a (linear) dipeptide was detected.
  7. Otherwise, check for cyclic dipeptides: search for a six-membered ring that contains two non‐aromatic peptide bonds.
  8. If none of these tests pass, return False.
  
If SMILES cannot be parsed, return (False, "Invalid SMILES string.").
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_dipeptide(smiles: str):
    """
    Determines if the given SMILES corresponds to a dipeptide.
    
    Our method accepts both linear dipeptides and cyclic dipeptides (diketopiperazines)
    by looking for valid peptide bonds (amide bonds that are not aromatic and not in rings)
    that, if "broken", yield two reasonably sized residues. For cyclic dipeptides we look 
    for a six-membered ring having exactly two such peptide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a dipeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # Check overall size: dipeptides should be relatively small.
    heavy_atoms_total = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(heavy_atoms_total) > 60:
        return False, f"Too many heavy atoms ({len(heavy_atoms_total)}) for a typical dipeptide."
    
    # Define SMARTS for basic amide (peptide) bond: carbonyl carbon single-bonded to nitrogen.
    amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")
    matches = mol.GetSubstructMatches(amide_pattern)
    if not matches:
        return False, "No amide (peptide) bond pattern detected."
    
    # Prepare to store candidate cleavage bonds.
    candidate_cleavages = []
    
    # Allowed heavy atom count per fragment for a typical amino acid residue.
    lower_bound = 4
    upper_bound = 35
    
    # Evaluate each candidate amide bond.
    for match in matches:
        # According to our SMARTS, first atom is carbonyl carbon, last is amide nitrogen.
        c_idx = match[0]
        n_idx = match[-1]
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond is None:
            continue
        # Exclude bonds if either atom is aromatic.
        if mol.GetAtomWithIdx(c_idx).GetIsAromatic() or mol.GetAtomWithIdx(n_idx).GetIsAromatic():
            continue
        # Exclude bonds that are in a ring (handle cyclic separately).
        if bond.IsInRing():
            continue
        
        # Attempt to cleave this bond.
        try:
            # FragmentOnBonds returns a molecule with dummy atoms ([*]) added at the cleavage sites.
            frag_mol = rdmolops.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        except Exception:
            continue
        
        frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        if len(frags) != 2:
            continue
        
        # Count heavy atoms (atomic number > 1) in each fragment.
        counts = []
        for frag in frags:
            count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
            counts.append(count)
        # Both fragments should have a size typical of residues.
        if all(lower_bound <= cnt <= upper_bound for cnt in counts):
            diff = abs(counts[0] - counts[1])
            candidate_cleavages.append((bond.GetIdx(), counts, diff))
    
    # If we found at least one candidate, choose the one with the smallest difference (most balanced)
    if candidate_cleavages:
        best = min(candidate_cleavages, key=lambda x: x[2])
        return True, ("Dipeptide detected (linear): selected peptide bond cleavage yields fragments "
                      f"with heavy atom counts {best[1][0]} and {best[1][1]}.")
    
    # Next, check for cyclic dipeptides: a diketopiperazine has a 6-membered ring containing 2 peptide bonds.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        bond_count_in_ring = 0
        for match in matches:
            c_idx = match[0]
            n_idx = match[-1]
            bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
            if bond is None:
                continue
            # Check that both atoms in the matched bond are in the ring.
            if c_idx in ring and n_idx in ring:
                # Exclude if aromatic.
                if mol.GetAtomWithIdx(c_idx).GetIsAromatic() or mol.GetAtomWithIdx(n_idx).GetIsAromatic():
                    continue
                bond_count_in_ring += 1
        if bond_count_in_ring == 2:
            return True, ("Cyclic dipeptide (diketopiperazine) detected: a 6-membered ring contains two peptide bonds.")
    
    return False, "No valid peptide bond connecting two appropriate fragments was found."

# Example usage with some test SMILES:
if __name__ == "__main__":
    test_smiles_list = [
        # True positives (dipeptides)
        "S(=O)(=N)(CC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)O",   # gamma-glutamylmethionine sulfoximine
        "SCC(NC(=O)C(N)C(O)C)C(O)=O",                         # Threoninyl-Cysteine
        "CC(C)[C@H](NC(=O)[C@H](C)N)C(O)=O",                   # Ala-Val
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCC(O)=O",            # Trp-Gly
        "C1C[C@@]2([H])C(N[C@H](C(N2C1)=O)CC(C)C)=O",         # cyclo(L-Leu-L-Pro) cyclic dipeptide
        # True negatives / ambiguous cases:
        "NCC(N1[C@@H](CCC1)C(=O)NC=2C=CC3=C(C2)C=CC=C3)=O",   # glycyl-L-proline 2-naphthylamide (ambiguous extra amide bonds)
    ]
    for sm in test_smiles_list:
        res, reason = is_dipeptide(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")