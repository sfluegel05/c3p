"""
Classifies: CHEBI:48030 tetrapeptide
"""
#!/usr/bin/env python3
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino‐acid residues connected by peptide linkages.

Improved Heuristic:
  - First, attempt to find a contiguous chain of four candidate alpha–carbons.
    We define an “alpha–carbon candidate” as a chiral carbon bound to an amino (N) and a carbonyl 
    (C(=O)), using two SMARTS patterns (for both possible chiralities).
  - In a linear peptide the expected graph distance (number of bonds) between consecutive 
    alpha–carbons is about 3.
  - As a fallback, we count peptide bonds using the backbone SMARTS ("C(=O)N[C@H]" and "C(=O)N[C@@H]").
  - In addition, overall properties (molecular weight and nitrogen count) are checked.
Note:
  This heuristic is still imperfect and might mis‐classify modified or cyclic peptides, or molecules 
  that accidentally contain a similar motif.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops
import itertools

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    
    The algorithm uses two complementary criteria:
      1. Try to identify a contiguous chain of 4 alpha–carbon candidates. 
         We define a candidate using a relaxed version of the alpha–residue motif:
         "[C@H](N)C(=O)" (and its mirror image). For a linear tetrapeptide, the graph distance
         (number of bonds) between consecutive alpha–carbons is expected to be 3.
      2. Otherwise, count peptide bonds by matching SMARTS for backbone fragments
         "C(=O)N[C@H]" and "C(=O)N[C@@H]". Though modified peptides may have extra amide bonds,
         a linear tetrapeptide typically has 3 (or sometimes 4) such bonds.
      
    Additional property checks include molecular weight and overall nitrogen atom count.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a tetrapeptide, False otherwise.
        str: Explanation for the decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pre-calculate the bond distance matrix (graph distances between atoms)
    dm = rdmolops.GetDistanceMatrix(mol)
    
    #############################
    # Criterion 1: Identify a contiguous chain of 4 alpha–carbons.
    # We use a relaxed SMARTS that requires a chiral carbon attached to an amino (N) group and a carbonyl.
    alpha_pat1 = Chem.MolFromSmarts("[C@H](N)C(=O)")
    alpha_pat2 = Chem.MolFromSmarts("[C@@H](N)C(=O)")
    alpha_candidates = set()
    if alpha_pat1:
        for match in mol.GetSubstructMatches(alpha_pat1):
            # match[0] is the chiral center (alpha–carbon)
            alpha_candidates.add(match[0])
    if alpha_pat2:
        for match in mol.GetSubstructMatches(alpha_pat2):
            alpha_candidates.add(match[0])
    
    # Try to find a contiguous chain of 4 alpha–carbons in which the graph distance between 
    # consecutive alpha–carbons is exactly 3 (a typical span in an amino acid backbone).
    found_chain = False
    chain = None
    alpha_candidates = list(alpha_candidates)
    if len(alpha_candidates) >= 4:
        # Check all permutations of 4 candidate atoms to find a contiguous chain.
        for combo in itertools.permutations(alpha_candidates, 4):
            if all(dm[combo[i], combo[i+1]] == 3 for i in range(3)):
                found_chain = True
                chain = combo
                break

    #############################
    # Criterion 2: Count peptide bonds using strict backbone fragment SMARTS.
    pat1 = Chem.MolFromSmarts("C(=O)N[C@H]")
    pat2 = Chem.MolFromSmarts("C(=O)N[C@@H]")
    peptide_bonds = set()
    if pat1:
        for match in mol.GetSubstructMatches(pat1):
            if len(match) >= 3:
                # Use tuple (carbonyl carbon idx, amide nitrogen idx)
                peptide_bonds.add((match[0], match[2]))
    if pat2:
        for match in mol.GetSubstructMatches(pat2):
            if len(match) >= 3:
                peptide_bonds.add((match[0], match[2]))
    n_peptide = len(peptide_bonds)
    
    #############################
    # Additional property checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Broadened molecular weight range to accommodate certain modified tetrapeptides.
    if not (250 <= mol_wt <= 2000):
        return False, "Molecular weight {:.1f} Da out of tetrapeptide range.".format(mol_wt)
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if not (4 <= n_nitrogen <= 30):
        return False, "Nitrogen count {} not typical for a tetrapeptide.".format(n_nitrogen)
    
    #############################
    # Combine the criteria.
    # First, if a contiguous chain of 4 alpha–carbons is found, we accept.
    if found_chain:
        return True, ("Found a contiguous chain of 4 alpha–carbons (indices: {}) and overall properties "
                       "consistent with a tetrapeptide.".format(chain))
    # Otherwise, fallback to peptide–bond counting.
    # Accept if we find either exactly 3 or 4 peptide bond fragments.
    elif n_peptide in [3, 4]:
        return True, ("Found {} peptide bond(s) matching the backbone motif and overall properties "
                       "consistent with a tetrapeptide.".format(n_peptide))
    else:
        return False, ("Found {} peptide bond(s) and {} alpha-residue candidate(s); expected either a contiguous "
                       "chain of 4 alpha–carbons or 3–4 peptide bonds for a linear tetrapeptide."
                       .format(n_peptide, len(alpha_candidates)))

# For quick testing, uncomment one or more test cases below:
if __name__ == '__main__':
    # Should be classified as tetrapeptide (true positive)
    tests = [
        ("Lys-Thr-Trp-Tyr",
         "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"),
        ("Asp-Phe-Val-Tyr",
         "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"),
        # A modified tetrapeptide that previously was a false negative:
        ("Z-DEVD-FMK",
         "N([C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(CF)=O)CC(=O)OC)=O)C(C)C)=O)CCC(OC)=O)=O)CC(OC)=O)C(OCC=1C=CC=CC1)=O"),
        # One potential false positive example:
        ("Pteropterin (expected non-tetrapeptide)",
         "C1=CC(=CC=C1C(=O)NC(CCC(=O)NC(CCC(=O)NC(CCC(=O)O)C(=O)O)C(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N")
    ]
    
    for name, smi in tests:
        res, reason = is_tetrapeptide(smi)
        print(f"{name}:\n  Result: {res}\n  Reason: {reason}\n")