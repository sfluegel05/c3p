"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine-containing compounds.
Definition: Compounds having the structure RC(=NR)NR2 (e.g. -C(=NH)NH2 and substituted variants).
This algorithm uses SMARTS for both neutral and cationic variants,
applies filters to exclude amidoxime-like moieties, verifies heavy-atom connectivity,
ensures proper sp2 hybridization for the carboxamidine carbon and double-bonded nitrogen,
and flags peptide-like molecules.
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
    
    # Define SMARTS patterns for the carboxamidine motif:
    # Pattern for neutral carboxamidine: carbon double-bonded to an [NX2] and single-bonded to an [NX3]
    pattern_neutral = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    # Pattern for a cationic variant (sometimes encountered)
    pattern_cationic = Chem.MolFromSmarts("[CX3](=[NX2+])[NX3]")
    if pattern_neutral is None or pattern_cationic is None:
        return False, "Error defining SMARTS patterns."
    
    # Gather all matches from both patterns.
    matches_neutral = mol.GetSubstructMatches(pattern_neutral)
    matches_cationic = mol.GetSubstructMatches(pattern_cationic)
    all_matches = set(matches_neutral) | set(matches_cationic)
    if not all_matches:
        return False, "Carboxamidine moiety (-C(=NH)NH2 or substituted variant) not found."
    
    valid_matches = []
    # Loop over each match.
    for match in all_matches:
        # Expect the order according to our SMARTS:
        #   match[0] = carboxamidine carbon,
        #   match[1] = double-bonded nitrogen,
        #   match[2] = single-bonded nitrogen.
        if len(match) < 3:
            continue
        idx_c, idx_ndbl, idx_nsingle = match[:3]
        atom_c = mol.GetAtomWithIdx(idx_c)
        atom_ndbl = mol.GetAtomWithIdx(idx_ndbl)
        atom_nsingle = mol.GetAtomWithIdx(idx_nsingle)
        
        # Filter A: Exclude if the double-bonded nitrogen is attached (by a single bond)
        # to any oxygen that carries at least one hydrogen, to avoid amidoxime groups.
        skip_match = False
        for nbr in atom_ndbl.GetNeighbors():
            # Skip atoms that are part of this motif.
            if nbr.GetIdx() in (idx_c, idx_nsingle):
                continue
            if nbr.GetAtomicNum() == 8:  # Oxygen
                bond = mol.GetBondBetweenAtoms(atom_ndbl.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                    if nbr.GetTotalNumHs() > 0:
                        skip_match = True
                        break
        if skip_match:
            continue
        
        # Filter B: Check that the carboxamidine carbon has the expected heavy-atom environment.
        # Count heavy (non-hydrogen) neighbors.
        heavy_nbrs = [nbr for nbr in atom_c.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_nbrs) not in (2, 3):
            continue
        # Exactly 2 of these must be nitrogen atoms.
        num_nitrogen = sum(1 for nbr in heavy_nbrs if nbr.GetAtomicNum() == 7)
        if num_nitrogen != 2:
            continue
        
        # Filter C: Ensure that the carboxamidine carbon is sp2-hybridized (for proper amidine geometry).
        if atom_c.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
            continue
        
        # Filter D: Ensure that the double-bonded nitrogen is sp2-hybridized.
        if atom_ndbl.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
            continue
        
        # If we reached here, we count this match as valid.
        valid_matches.append(match)
    
    if not valid_matches:
        return False, "Carboxamidine moiety not found after refined filtering."
    
    # Optionally, reject molecules that appear peptide-like:
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if num_amide >= 4 and mol_wt > 600:
        return False, "Molecule appears peptide-like (many amide bonds in a heavy molecule), likely a false positive."
    
    return True, f"Found carboxamidine group in {len(valid_matches)} location(s)."

# For testing/debugging (run only when this script is executed directly)
if __name__ == "__main__":
    test_smiles = [
        "[H]C(N)=N",  # formamidine: expected true positive.
        "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",  # (Z)-acetamiprid.
        "CC(=N)NCC1=CC(CN)=CC=C1",  # N-[3-(aminomethyl)benzyl]acetamidine.
        "CN1CCCN=C1\\C=C\\c1cccs1",  # pyrantel.
        "CC1=N[C@@H]([C@@H](O)CN1)C(O)=O",  # 5-hydroxyectoine.
        "CN(Cc1ccc(Cl)nc1)C(C)=NC#N",  # acetamiprid.
        "CCCCCCOC(=O)\\N=C(\\N)c1ccc(NCc2nc3cc(ccc3n2C)C(=O)N(CCC(=O)OCC)c2ccccn2)cc1",  # dabigatran etexilate.
        "[H][C@]1(CCN1C(=O)[C@H](NCC(O)=O)C1CCCCC1)C(=O)NCc1ccc(cc1)C(N)=N",  # melagatran.
        "NC(=N)c1ccc(NN=Nc2ccc(cc2)C(N)=N)cc1",  # diminazene.
        "NC(N)=N",  # guanidine (expected to be classified as carboxamidine).
        # A known false positive from previous attempt:
        "N(CC1=CC=CC=C1)(C/C(=N/O)/N)C2=CC=CC=C2",  # cetoxime.
        # A peptide-like false positive:
        "O=C(N[C@@H](CC(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)CCCN=C(N)N",  # Val-Arg-Leu peptide fragment.
        # Previously missed cases:
        "NC(=N)N1CC(O)c2ccccc2C1",  # 4-hydroxydebrisoquin.
        "NC(=N)N1CCc2ccccc2C1",  # debrisoquin.
        "CCCCCNC(=N)N\\N=C\\c1c[nH]c2ccc(OC)cc12",  # tegaserod.
        "CNC(NCCSCc1csc(CN(C)C)n1)=C[N+]([O-])=O",  # nizatidine.
    ]
    
    for smi in test_smiles:
        res, reason = is_carboxamidine(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")