"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: Oligopeptide – a peptide containing a relatively small number of amino acids.

This version uses modified heuristics:
  1. It searches for amide (peptide) bonds using the SMARTS "C(=O)N" and then filters out those
     that are in rings (often seen in unrelated chemistry).
  2. It estimates the residue count as (# peptide bonds + 1) and accepts between 2 and 10 residues.
  3. It counts chiral alpha–carbon units using the SMARTS "[C@H](N)" and "[C@@H](N)".
  4. It checks that all peptide-bond atoms (carbonyl carbon and the adjacent amide nitrogen)
     lie in a single fragment.
  5. For 3 or more residues it requires that a contiguous backbone appears
     (via substructure searches for "N[C@H](*)C(=O)N" or its mirror).
  6. It enforces a molecular weight range roughly estimated as (n_residues * 60) < MW < (n_residues * 200)
     (this may help rule out compounds with peptide‐like bonds that are too heavy).
     
Note: This is a heuristic; many molecules are borderline and some failures (false positives/negatives)
may still occur.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines whether the given molecule (via its SMILES) is a small peptide (oligopeptide).

    Heuristics used:
      - Search for amide (peptide) bonds (SMARTS "C(=O)N"), but ignore those that are in rings.
      - Estimate residue count as (number of peptide bonds + 1); allow between 2 and 10 residues.
      - Count alpha–carbon centers (matches to "[C@H](N)" and "[C@@H](N)").
         For dipeptides, require exactly 2; for larger chains, require at least 2.
      - Ensure all peptide bond atoms (carbonyl C and amide N) lie in a single fragment.
      - For peptides with ≥3 residues, require that a backbone motif is found, for example
            "N[C@H](*)C(=O)N" or "N[C@@H](*)C(=O)N".
      - Check that the molecular weight roughly falls in the expected range:
            lower bound = n_residues * 60 Da, upper bound = n_residues * 200 Da.
      
    Args:
        smiles (str): SMILES string of the molecule.
      
    Returns:
        bool: True if the molecule is classified as an oligopeptide, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Find "peptide" (amide) bonds via SMARTS "C(=O)N"
    pb_smarts = Chem.MolFromSmarts("C(=O)N")
    all_pb_matches = mol.GetSubstructMatches(pb_smarts)
    
    # Only consider peptide bonds whose carbonyl C and the adjacent N are not in a ring.
    pb_matches = []
    for match in all_pb_matches:
        # In our SMARTS, match[0] is the carbonyl carbon; match[2] is the attached N.
        # (match[1] corresponds to the oxygen, but we ignore it here)
        c_atom = mol.GetAtomWithIdx(match[0])
        n_atom = mol.GetAtomWithIdx(match[2])
        if c_atom.IsInRing() or n_atom.IsInRing():
            continue
        pb_matches.append(match)
    
    n_pbonds = len(pb_matches)
    if n_pbonds == 0:
        return False, "No peptide (amide) bonds found (or none outside ring systems)"
    
    n_residues = n_pbonds + 1
    if n_residues < 2:
        return False, "Too few peptide bonds to form a peptide (need at least 2 residues)"
    if n_residues > 10:
        return False, f"Found {n_residues} amino acid residues which is too many for an oligopeptide"

    # Step 2. Count alpha–carbon centers (with declared chirality) using SMARTS.
    alpha1 = Chem.MolFromSmarts("[C@H](N)")
    alpha2 = Chem.MolFromSmarts("[C@@H](N)")
    alpha_set = set()
    for m in mol.GetSubstructMatches(alpha1):
        # m[0] is the chiral carbon.
        alpha_set.add(m[0])
    for m in mol.GetSubstructMatches(alpha2):
        alpha_set.add(m[0])
    n_alpha = len(alpha_set)
    if n_residues == 2 and n_alpha != 2:
        return False, f"For a dipeptide, exactly 2 alpha–carbon centers are expected; found {n_alpha}"
    if n_residues >= 3 and n_alpha < 2:
        return False, f"Expected at least 2 alpha–carbon centers; found {n_alpha}"
    
    # Step 3. Check that all atoms involved in peptide bonds are in one contiguous fragment.
    pb_atom_idxs = set()
    for m in pb_matches:
        pb_atom_idxs.add(m[0])  # carbonyl carbon
        pb_atom_idxs.add(m[2])  # amide nitrogen
    frags = Chem.GetMolFrags(mol, asMols=False)
    count_frag_with_pb = sum(1 for frag in frags if any(a in frag for a in pb_atom_idxs))
    if count_frag_with_pb > 1:
        return False, "Peptide-bond atoms are not contained in a single contiguous fragment"
    
    # Step 4. For peptides with ≥3 residues, look for a contiguous backbone motif.
    if n_residues >= 3:
        bb_smarts1 = Chem.MolFromSmarts("N[C@H](*)C(=O)N")
        bb_smarts2 = Chem.MolFromSmarts("N[C@@H](*)C(=O)N")
        if not (mol.HasSubstructMatch(bb_smarts1) or mol.HasSubstructMatch(bb_smarts2)):
            return False, "Peptide-bond connectivity does not appear contiguous (backbone motif missing)"
    
    # Step 5. Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    min_expected = n_residues * 60
    max_expected = n_residues * 200  # a bit looser than before to allow variability
    if mol_wt < min_expected:
        return False, (f"Molecular weight ({mol_wt:.1f} Da) is too low for a {n_residues}-residue peptide "
                       f"(expected at least {min_expected} Da)")
    if mol_wt > max_expected:
        return False, (f"Molecular weight ({mol_wt:.1f} Da) is too high for a {n_residues}-residue peptide "
                       f"(expected at most {max_expected} Da)")
    
    # Additional property: count rotatable bonds.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    reason = (f"Detected {n_pbonds} peptide bond(s) (≈{n_residues} residue(s)), {n_rotatable} rotatable bond(s), "
              f"{n_alpha} alpha–carbon center(s), and MW of {mol_wt:.1f} Da. "
              "Peptide bond atoms are in a single fragment and a contiguous backbone motif is detected. "
              "This is consistent with an oligopeptide.")
    
    return True, reason

# For local testing one might do:
if __name__ == "__main__":
    # Example: Leu-Trp dipeptide
    test_smiles = "CC(C)C[C@H](N)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O"
    result, explanation = is_oligopeptide(test_smiles)
    print(result, explanation)