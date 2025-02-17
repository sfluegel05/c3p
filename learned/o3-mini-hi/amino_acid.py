"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies chemical entities of the class amino acid.
Definition used: An amino acid is a molecule containing at least one carboxylic 
acid group (in neutral or deprotonated form) and one amino group that is not 
acylated (“free”). Moreover, at least one acid–amine pair must occur “close together”
(i.e. with a shortest bond path of 4 bonds or fewer) and without any intervening peptide (amide) bonds,
so that the two functional groups belong to the same amino acid “residue.”
Note: Some borderline cases may still be ambiguous.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is (likely) an amino acid based on its SMILES string.
    The algorithm requires:
      1. Valid parsing.
      2. At least one carboxylic acid group (neutral –C(=O)[OH] or anionic –C(=O)[O-]).
      3. At least one free (non-acylated) amino group.
      4. The existence of at least one acid–amine pair with a shortest bond path of 4 bonds or fewer.
         Additionally, we inspect the bond–path: if any bond in the path is part of an amide bond
         (as seen in di- or oligopeptides) then this pair is rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an amino acid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- 1. Identify carboxyl groups ---
    # We use two SMARTS patterns to capture the acid in its neutral and anionic forms.
    acid_neutral_smarts = "[CX3](=O)[OX2H]"
    acid_anion_smarts   = "[CX3](=O)[O-]"
    acid_neutral = Chem.MolFromSmarts(acid_neutral_smarts)
    acid_anion   = Chem.MolFromSmarts(acid_anion_smarts)
    
    acid_indices = set()
    for match in mol.GetSubstructMatches(acid_neutral):
        # match[0] is the carbon of the acid group
        acid_indices.add(match[0])
    for match in mol.GetSubstructMatches(acid_anion):
        acid_indices.add(match[0])
        
    if not acid_indices:
        return False, "No carboxylic acid group found"

    # --- 2. Identify free (non-acylated) amino groups ---
    # Instead of using a fixed SMARTS we iterate over all nitrogen atoms and discard those directly
    # bonded to a carbonyl carbon (except when the carbonyl is part of a proper acid group).
    free_amino_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        n_idx = atom.GetIdx()
        is_acylated = False
        # Check every neighbor: if the nitrogen is directly bonded to a carbon that is part of an amide linkage,
        # then mark it as acylated.
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() != "C":
                continue
            c_idx = neighbor.GetIdx()
            # Examine bonds from this carbon to its neighbors (other than our N)
            # Count bonds to O in double bond; if the C also has a single-bonded O with H or - charge, it is likely
            # a carboxyl group (acid) rather than an amide.
            double_o = False
            o_single = False
            for nbr in neighbor.GetNeighbors():
                if nbr.GetIdx() == n_idx:
                    continue
                if nbr.GetSymbol() == "O":
                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx())
                    if bond is None:
                        continue
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_o = True
                    elif bond.GetBondType() == Chem.BondType.SINGLE:
                        # Check if this oxygen is part of an acid (has a hydrogen or negative charge)
                        # (Note: using formal charge here; if O has +1 charge it’s not acid)
                        if nbr.GetTotalNumHs() >= 1 or nbr.GetFormalCharge() < 0:
                            o_single = True
            # If the C has a double-bonded O but NOT also a single-bonded oxygen typical for a carboxyl,
            # we take it as part of an amide bond.
            if double_o and not o_single:
                is_acylated = True
                break
        if not is_acylated:
            free_amino_indices.append(n_idx)
    
    if not free_amino_indices:
        return False, "No free (non-acylated) amino group found"
        
    # --- 3. Exclude molecules that are likely peptides ---
    # Many di- or oligopeptides will also have one free acid and one free amino group.
    # (For our purposes, we only want amino acid “residues”.)
    # Here we use the heuristic that if the free acid and free amine do not belong to the same fragment,
    # then we reject.
    # Also, if none of the acid–amine pairs appears to come from the same amino acid residue, we reject.
    for acid_idx in acid_indices:
        acid_fragment = rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False)[0]
        # (In RDKit GetMolFrags returns a tuple of tuples of atom indices; if the acid atom appears
        # in one fragment, then the candidate amino group must lie in the same fragment.)
        for amine_idx in free_amino_indices:
            # Ensure that the acid and amino fall in the same connected component.
            frags = rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
            frag_of_acid = None
            frag_of_amine = None
            for frag in frags:
                if acid_idx in frag:
                    frag_of_acid = frag
                if amine_idx in frag:
                    frag_of_amine = frag
            if frag_of_acid is None or frag_of_amine is None or frag_of_acid != frag_of_amine:
                continue  # not in the same fragment
            # Compute the shortest path (list of atom indices) between acid carbon and amine nitrogen.
            path = rdmolops.GetShortestPath(mol, acid_idx, amine_idx)
            if not path or len(path) < 2:
                continue
            # The bond distance is len(path)-1.
            distance = len(path) - 1
            # We require a short distance (≤4 bonds) as a necessary heuristic.
            if distance > 4:
                continue

            # --- 4. Inspect intervening bonds to see if a peptide (amide) bond is present ---
            # For each bond along the path, when we see a bond between a carbon and nitrogen,
            # we check if that C bears a double-bonded oxygen that is not paired with an -OH (which in acids is expected).
            peptide_bond_found = False
            # Loop over consecutive atom pairs in the path.
            for i in range(len(path)-1):
                a1 = mol.GetAtomWithIdx(path[i])
                a2 = mol.GetAtomWithIdx(path[i+1])
                # Look only at C–N bonds (in either order)
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "N"}:
                    # Identify the C atom.
                    c_atom = a1 if a1.GetSymbol() == "C" else a2
                    # Check bonds from c_atom to oxygen.
                    double_o = False
                    single_o = False
                    for nbr in c_atom.GetNeighbors():
                        if nbr.GetSymbol() == "O":
                            bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                            if bond is None:
                                continue
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                double_o = True
                            elif bond.GetBondType() == Chem.BondType.SINGLE:
                                if nbr.GetTotalNumHs() >= 1 or nbr.GetFormalCharge() < 0:
                                    single_o = True
                    # If c_atom is bonded to a double-bonded oxygen but not the pattern expected of a carboxyl group (i.e. lacking an -OH),
                    # we flag this bond as amide-like.
                    if double_o and not single_o:
                        peptide_bond_found = True
                        break
            if peptide_bond_found:
                # This acid–amine pair is connected but through an intervening peptide-bond;
                # so we skip this pair.
                continue
            else:
                # We found a candidate acid–amine pair that are close (≤4 bonds) and do not include an amide bond.
                return True, f"Found free carboxylic acid and free amino group with a bond distance of {distance} (within cutoff)."
                
    # If no acid–amine pair meeting all the criteria is found then reject.
    return False, "Functional groups found but none are close enough (or belong to the same amino acid residue) without peptide linkages."

# Example test cases (can be removed in production code)
if __name__ == "__main__":
    test_cases = {
        "(S)-gabaculine": "C1=C(C[C@@H](C=C1)N)C(O)=O",
        "porphyra-334": "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",
        "L-thyroxine": "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",
        "robenacoxib": "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",
        # A dipeptide (false positive in previous version) — should be rejected:
        "Glu-Cys-Thr (peptide-like)": "SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)C(O)=O)CCC(=O)N)C(O)=O",
        # A false negative test: N,N-dihydroxydihomomethionine should now be accepted.
        "N,N-dihydroxydihomomethionine": "CSCCCCC(N(O)O)C(O)=O",
    }
    for name, sm in test_cases.items():
        result, reason = is_amino_acid(sm)
        print(f"{name}: {result} ({reason})")