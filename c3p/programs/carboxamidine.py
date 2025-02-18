"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine-containing compounds.
Definition: Compounds having the structure RC(=NR)NR2, i.e. a -C(=NH)NH2 group or a substituted variant.
This improved algorithm refines the SMARTS match and filtering criteria:
  • Uses both neutral and cationic SMARTS for the carboxamidine motif.
  • Excludes cases where the double‐bonded N is attached to an oxygen bearing at least one hydrogen.
  • Checks that the carboxamidine carbon (after adding hydrogens) has a proper heavy-atom environment:
       – Its heavy (non‐H) neighbors should total either 2 (if R is H) or 3 (if R is not hydrogen)
         and exactly two of these must be nitrogen atoms.
  • Rejects molecules that appear peptide-like (many amide bonds in a heavy molecule).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2).

    Args:
      smiles (str): SMILES string representing the molecule.

    Returns:
      bool: True if a carboxamidine group is detected, False otherwise.
      str: Explanation of the classification outcome.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    # Add explicit hydrogens so that neighbor analysis is unambiguous.
    mol = Chem.AddHs(mol)
    
    # Define two SMARTS patterns for carboxamidine:
    pattern_neutral = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    pattern_cationic = Chem.MolFromSmarts("[CX3](=[NX2+])[NX3]")
    if pattern_neutral is None or pattern_cationic is None:
        return False, "Error defining SMARTS patterns."
    
    # Get all substructure matches from both patterns.
    matches_neutral = mol.GetSubstructMatches(pattern_neutral)
    matches_cationic = mol.GetSubstructMatches(pattern_cationic)
    all_matches = set(matches_neutral) | set(matches_cationic)
    if not all_matches:
        return False, "Carboxamidine moiety (-C(=NH)NH2 or substituted equivalent) not found."
    
    valid_matches = []
    for match in all_matches:
        # According to our SMARTS match order we expect:
        #   match[0] = carboxamidine carbon,
        #   match[1] = double-bonded nitrogen,
        #   match[2] = single-bonded nitrogen.
        if len(match) < 3:
            continue
        idx_c, idx_ndbl, idx_nsingle = match[:3]
        atom_c = mol.GetAtomWithIdx(idx_c)
        atom_ndbl = mol.GetAtomWithIdx(idx_ndbl)
        atom_nsingle = mol.GetAtomWithIdx(idx_nsingle)
        
        # Filter A: Exclude if the double-bonded nitrogen is attached (via a single bond)
        # to any oxygen that carries at least one hydrogen (i.e. avoid amidoxime groups).
        skip_match = False
        for nbr in atom_ndbl.GetNeighbors():
            # Skip atoms which are part of the carboxamidine motif.
            if nbr.GetIdx() in (idx_c, idx_nsingle):
                continue
            if nbr.GetAtomicNum() == 8:  # Oxygen
                bond = mol.GetBondBetweenAtoms(idx_ndbl, nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    if nbr.GetTotalNumHs() > 0:
                        skip_match = True
                        break
        if skip_match:
            continue
        
        # Filter B: Check that the carboxamidine carbon has exactly the expected heavy-atom environment.
        # Get heavy (non-hydrogen) neighbors of the carbon.
        heavy_nbrs = [nbr for nbr in atom_c.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # In a true carboxamidine group, the carboxamidine carbon is bound to the two N atoms of the motif
        # plus one extra substituent. When R = H the heavy-neighbor count is 2 (only nitrogens),
        # when R ≠ H it is 3. So acceptable counts are 2 or 3 and exactly 2 must be N.
        if len(heavy_nbrs) not in (2, 3):
            continue
        num_nitrogen = sum(1 for nbr in heavy_nbrs if nbr.GetAtomicNum() == 7)
        if num_nitrogen != 2:
            continue
        
        # Filter C: Ensure the carboxamidine carbon is exocyclic (not aromatic).
        if atom_c.GetIsAromatic():
            continue

        # (Optional) Filter D: Ensure that the double-bonded N is sp2-hybridized.
        if atom_ndbl.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
            continue

        valid_matches.append(match)
    
    if not valid_matches:
        return False, "Carboxamidine moiety not found after refined filtering."
    
    # In heavy molecules, many amide bonds are a good indication of peptide-like systems.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if num_amide >= 4 and mol_wt > 600:
        return False, "Molecule appears to be peptide-like (many amide bonds in a heavy molecule), likely a false positive."
    
    return True, f"Found carboxamidine group in {len(valid_matches)} location(s)."

# For testing/debugging (run only when this script is executed directly)
if __name__ == "__main__":
    test_smiles = [
        "[H]C(N)=N",  # formamidine: expected true positive.
        "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",  # (Z)-acetamiprid.
        "CC(=N)NCC1=CC(CN)=CC=C1",  # N-[3-(aminomethyl)benzyl]acetamidine.
        "CN1CCCN=C1\\C=C\\c1cccs1",  # pyrantel.
        "CC1=N[C@@H]([C@@H](O)CN1)C(O)=O",  # 5-hydroxyectoine.
        "CNC(N)=N",  # acetamidine.
        "NC(=N)c1ccccc1",  # benzamidine.
        # A known false positive from previous attempt:
        "N(CC1=CC=CC=C1)(C/C(=N/O)/N)C2=CC=CC=C2",  # cetoxime.
        # A peptide-like false positive:
        "O=C(N[C@@H](CC(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)CCCN=C(N)N",  # Val-Arg-Leu.
        # A known false negative previously:
        "CNC(NCCSCc1csc(CN(C)C)n1)=C[N+]([O-])=O",  # nizatidine.
    ]
    
    for smi in test_smiles:
        res, reason = is_carboxamidine(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")