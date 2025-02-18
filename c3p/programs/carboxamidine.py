"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine-containing compounds.
Definition: Compounds having the structure RC(=NR)NR2, i.e. the –C(=NH)NH2 group or a substituted variant.
The algorithm improves on the previous version by refining the SMARTS match and filtering criteria:
  • It uses two SMARTS patterns (neutral and cationic) for the carboxamidine motif.
  • It excludes cases where the double‐bonded nitrogen is bound to an oxygen that carries at least one hydrogen (i.e. amidoxime moieties).
  • It checks that the carboxamidine carbon has exactly three neighbors (two nitrogens and one extra substituent) and is not aromatic,
    which is expected for an exocyclic carboxamidine group.
  • It uses an amide-bond count in heavy molecules (mol wt >600) to flag systems that are peptide-like.
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
    # Create molecule from SMILES and add explicit hydrogens for correct neighbor analysis
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # Define SMARTS patterns for the carboxamidine motif.
    # The neutral pattern: carbon (3-connected) double bonded to a nitrogen with two connections, 
    # and a single-bonded nitrogen.
    pattern_neutral = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    pattern_cationic = Chem.MolFromSmarts("[CX3](=[NX2+])[NX3]")
    if pattern_neutral is None or pattern_cationic is None:
        return False, "Error defining SMARTS patterns."
    
    # Get all substructure matches (as tuples of atom indices) from both patterns.
    matches1 = mol.GetSubstructMatches(pattern_neutral)
    matches2 = mol.GetSubstructMatches(pattern_cationic)
    all_matches = set(matches1) | set(matches2)
    if not all_matches:
        return False, "Carboxamidine moiety (-C(=NH)NH2 or substituted equivalent) not found."
    
    valid_matches = []
    for match in all_matches:
        # Expect match order: (idx_c, idx_dblN, idx_singleN)
        idx_c, idx_ndbl, idx_nsingle = match
        atom_c = mol.GetAtomWithIdx(idx_c)
        atom_ndbl = mol.GetAtomWithIdx(idx_ndbl)
        atom_nsingle = mol.GetAtomWithIdx(idx_nsingle)
        
        # Filter (A): Exclude if the double bonded N is attached (via a single bond) to any oxygen bearing at least one hydrogen.
        skip_match = False
        for nbr in atom_ndbl.GetNeighbors():
            # Do not check bonds to the carboxamidine carbon or the single bonded N
            if nbr.GetIdx() in (idx_c, idx_nsingle):
                continue
            if nbr.GetAtomicNum() == 8:  # Oxygen
                bond = mol.GetBondBetweenAtoms(idx_ndbl, nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    # if oxygen carries one or more explicit hydrogen, skip this match.
                    if nbr.GetTotalNumHs() > 0:
                        skip_match = True
                        break
        if skip_match:
            continue
        
        # Filter (B): Ensure the carboxamidine carbon has exactly three neighbors (it should be bound to two N atoms and a single R group).
        # (Our use of Chem.AddHs means that implicit hydrogens are now explicit.)
        if len(atom_c.GetNeighbors()) != 3:
            continue
        
        # Filter (C): The carboxamidine carbon is expected to be exocyclic so it should not be marked as aromatic.
        if atom_c.GetIsAromatic():
            continue
        
        valid_matches.append(match)
    
    if not valid_matches:
        return False, "Carboxamidine moiety not found after refined filtering."
    
    # As in the previous version: if the molecule is heavy and has many amide bonds it might be peptide-like.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if num_amide >= 4 and mol_wt > 600:
        return False, "Molecule appears to be peptide-like (many amide bonds in a heavy molecule), likely a false positive."
    
    return True, f"Found carboxamidine group in {len(valid_matches)} location(s)."

# For testing/debugging only (executed when running the script directly)
if __name__ == "__main__":
    # A few examples (many more were provided in the original outcomes)
    test_smiles = [
        "[H]C(N)=N",  # formamidine: expected true positive.
        "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",  # (Z)-acetamiprid
        "CC(=N)NCC1=CC(CN)=CC=C1",       # N-[3-(aminomethyl)benzyl]acetamidine
        "CN1CCCN=C1\\C=C\\c1cccs1",       # pyrantel
        "NC(=N)N1CC(O)c2ccccc2C1",        # 4-hydroxydebrisoquin (previous false negative)
        "CNC(N)=N",                     # guanidine (previous false negative)
        "NC(=N)c1ccccc1",               # benzamidine, expected true positive.
    ]
    
    for smi in test_smiles:
        res, reason = is_carboxamidine(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")