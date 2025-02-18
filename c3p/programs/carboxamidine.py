"""
Classifies: CHEBI:35359 carboxamidine
"""
#!/usr/bin/env python
"""
Classifies: Carboxamidine-containing compounds.
Definition: Compounds having the structure RC(=NR)NR2,
where the group –C(=NH)NH2 or a substituted analogue is present.
Examples include formamidine, acetamidine, benzamidine, guanidine etc.
The algorithm refines matching by:
  • Trying two SMARTS patterns (one for neutral and one for cationic representations)
  • Filtering out matches where the double‐bonded N is bound to oxygen (to avoid amidoximes)
  • Counting amide bonds and molecular weight to flag peptide-like false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2).
    
    To improve accuracy:
      - We try two SMARTS patterns to capture both neutral and positively charged forms.
      - For each match, we check that the double-bonded nitrogen is not bound to any oxygen.
      - Finally, we count amide bonds; if there are two or more amide bonds in a heavy molecule,
        we assume it is peptide-like and likely a false positive.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as containing a carboxamidine group.
        str: Explanation for the classification (or rejection).
    """
    # Parse input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # (For neighbor analysis, include implicit hydrogens as explicit)
    mol = Chem.AddHs(mol)
    
    # Define two SMARTS patterns:
    # pattern1: neutral representation: a trigonal carbon double-bonded to a (normally coordinated) nitrogen
    pattern1 = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    # pattern2: allow possibility that the double-bonded nitrogen carries a positive charge.
    pattern2 = Chem.MolFromSmarts("[CX3](=[NX2+])[NX3]")
    if pattern1 is None or pattern2 is None:
        return False, "Error defining SMARTS patterns."
    
    # Get all matches from both patterns (using sets to avoid duplicates).
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    all_matches = set(matches1) | set(matches2)
    
    if not all_matches:
        return False, "Carboxamidine moiety (-C(=NH)NH2 or substituted equivalent) not found."
    
    valid_matches = []
    # For each match, the SMARTS places atoms as:
    # index0: central carbon; index1: double-bonded N; index2: single-bonded N.
    for match in all_matches:
        idx_c, idx_ndbl, idx_nsingle = match
        atom_ndbl = mol.GetAtomWithIdx(idx_ndbl)
        # Filter out match if the double-bonded N is bound to any oxygen (e.g. amidoxime)
        # We check all neighbors of the double-bonded nitrogen.
        has_oxygen = any(neigh.GetAtomicNum() == 8 for neigh in atom_ndbl.GetNeighbors())
        if has_oxygen:
            continue
        valid_matches.append(match)
    
    if not valid_matches:
        return False, "Carboxamidine moiety not found after filtering out oxygen-substituted forms."
    
    # Count amide bonds in the molecule. This helps flag peptides that contain multiple C(=O)N bonds.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # If the molecule is large and contains several amide bonds, it is likely a peptide.
    if num_amide >= 2 and mol_wt > 300:
        return False, "Molecule appears to be peptide-like (many amide bonds in a heavy molecule), likely a false positive."
    
    return True, f"Found carboxamidine group in {len(valid_matches)} location(s)."

# For testing/debugging purposes (run only when executed as a script):
if __name__ == "__main__":
    # List of example SMILES (both true positives and some known problematic ones)
    test_smiles = [
        "[H]C(N)=N",  # formamidine (true positive)
        "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",  # (Z)-acetamiprid (true positive)
        "CC(=N)NCC1=CC(CN)=CC=C1",  # N-[3-(aminomethyl)benzyl]acetamidine (true positive)
        "CN1CCCN=C1\\C=C\\c1cccs1",  # pyrantel (true positive)
        "CC1=N[C@@H]([C@@H](O)CN1)C(O)=O",  # 5-hydroxyectoine (true positive)
        "CN(Cc1ccc(Cl)nc1)C(C)=NC#N",  # acetamiprid (true positive)
        "CCCCCCOC(=O)\\N=C(\\N)c1ccc(NCc2nc3cc(ccc3n2C)C(=O)N(CCC(=O)OCC)c2ccccn2)cc1",  # dabigatran etexilate
        # False positive examples (peptide-like):
        "N(CC1=CC=CC=C1)(C/C(=N/O)/N)C2=CC=CC=C2",  # cetoxime -- should be rejected
        "O=C(N[C@@H](CCCN=C(N)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)CCCN=C(N)N",  # Val-Arg-Leu peptide fragment
        # False negative example:
        "CNC(NCCSCc1csc(CN(C)C)n1)=C[N+]([O-])=O"  # nizatidine - should be accepted now
    ]
    
    for smi in test_smiles:
        res, reason = is_carboxamidine(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")