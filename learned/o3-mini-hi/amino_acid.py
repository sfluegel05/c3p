"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies chemical entities of the class amino acid.
Definition used: An amino acid is a molecule containing at least one free (non‐acylated) amino group and at least one carboxylic acid group 
(in neutral –C(=O)[OH] or anionic –C(=O)[O-] form). In addition, at least one acid–amine pair must appear “close together” 
(i.e. with a shortest bond path of 4 bonds or fewer) – ideally sharing an “alpha‐carbon” (i.e. acid carbon connected directly 
to a carbon that is also bonded to the free amino nitrogen).
Some borderline cases may still be ambiguous.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is (likely) an amino acid based on its SMILES string.

    The algorithm is as follows:
      1. Parse the molecule.
      2. Find all carboxylic acid groups (either neutral –C(=O)[OH] or anionic –C(=O)[O-]) where the acid carbon has degree=3.
      3. Find all “free” amino groups (nitrogen atoms that (a) are not acylated –i.e. not directly bonded to a carbonyl that lacks –OH,
         and (b) have zero formal charge). (We relax the condition on number of heavy neighbors to allow slightly substituted amines.)
      4. First attempt: Look for an alpha‐amino acid core defined as a carboxyl carbon bonded to an sp³ carbon (the candidate “α‐carbon”)
         that in turn is bonded to a free amino nitrogen. If found, return success with bond path length = 2.
      5. If not found, then for each acid–amine pair that lie in the same connected fragment, compute the shortest bond path.
         Accept the pair only if the bond distance is ≤ 4 and (if the distance is >2) at least one intervening atom (candidate “α‐carbon”) is sp³.
         Also verify that no bond along the path looks “amide‐like” (i.e. a C–N bond where the C is carbonyl and not –OH substituted).
    
    Args:
       smiles (str): SMILES string of the molecule.
        
    Returns:
       bool: True if the molecule is classified as an amino acid, False otherwise.
       str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # -------- Step 1: Find acid groups --------
    acid_neutral_smarts = "[CX3](=O)[OX2H]"
    acid_anion_smarts   = "[CX3](=O)[O-]"
    acid_neutral = Chem.MolFromSmarts(acid_neutral_smarts)
    acid_anion   = Chem.MolFromSmarts(acid_anion_smarts)
    acid_indices = set()  # store indices of acid carbons (the carbonyl in the acid group)

    # Look for neutral acids
    for match in mol.GetSubstructMatches(acid_neutral):
        # match: (acid carbon, hydroxyl oxygen) – take acid carbon:
        acid_atom = mol.GetAtomWithIdx(match[0])
        if acid_atom.GetDegree() == 3:
            acid_indices.add(match[0])
    # And anionic acids
    for match in mol.GetSubstructMatches(acid_anion):
        acid_atom = mol.GetAtomWithIdx(match[0])
        if acid_atom.GetDegree() == 3:
            acid_indices.add(match[0])
    if not acid_indices:
        return False, "No carboxylic acid group found"

    # -------- Step 2: Identify free (non‑acylated) amino groups --------
    free_amino_indices = []  # list of indices of N atoms that are free
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        if atom.GetFormalCharge() != 0:
            continue

        # Check if this nitrogen is acylated.
        # If any neighbor that is carbon is bonded to a double-bonded oxygen (and not to an –OH), then assume acylation.
        is_acylated = False
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            has_carbonyl = False
            has_hydroxyl = False
            for nn in nbr.GetNeighbors():
                if nn.GetIdx() == atom.GetIdx():
                    continue
                if nn.GetSymbol() == "O":
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                    if bond is None:
                        continue
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        has_carbonyl = True
                    elif bond.GetBondType() == Chem.BondType.SINGLE:
                        # if oxygen has at least one hydrogen (or a negative charge), it may be –OH.
                        if nn.GetTotalNumHs() >= 1 or nn.GetFormalCharge() < 0:
                            has_hydroxyl = True
            if has_carbonyl and not has_hydroxyl:
                is_acylated = True
                break
        if not is_acylated:
            free_amino_indices.append(atom.GetIdx())

    if not free_amino_indices:
        return False, "No free (non-acylated) amino group found"

    # -------- Step 3: First try to find an alpha-amino acid motif via direct connectivity --------
    # For each acid group, see if the acid carbon is bonded to a carbon (candidate alpha-carbon)
    # that in turn is bonded to a free amine nitrogen.
    for acid_idx in acid_indices:
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            # candidate alpha-carbon must preferably be sp3
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Check if any neighbor of candidate (other than the acid carbon) is a free amine.
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == acid_idx:
                    continue
                if nbr2.GetSymbol() == "N" and nbr2.GetIdx() in free_amino_indices:
                    return True, "Found acid carbon bonded to an sp3 alpha-carbon which is directly bonded to a free amino group (bond distance = 2)."

    # -------- Step 4: Relaxed search among acid and free amine groups via shortest bond paths --------
    # Get connected fragments to ensure acid and amine are in the same part of the molecule.
    frags = list(rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False))
    for acid_idx in acid_indices:
        for amine_idx in free_amino_indices:
            # Verify acid and amine are in same fragment.
            acid_frag = None
            amine_frag = None
            for frag in frags:
                if acid_idx in frag:
                    acid_frag = frag
                if amine_idx in frag:
                    amine_frag = frag
            if acid_frag is None or amine_frag is None or acid_frag != amine_frag:
                continue

            # Compute shortest path (list of atom indices).
            path = rdmolops.GetShortestPath(mol, acid_idx, amine_idx)
            if not path or len(path) < 2:
                continue
            distance = len(path) - 1  # number of bonds
            if distance > 4:
                continue  # too far apart to be part of one amino acid residue

            # If distance > 2, ensure that at least one of the intermediate atoms is sp3 (candidate alpha-carbon).
            if distance > 2:
                found_sp3 = False
                for idx in path[1:-1]:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        found_sp3 = True
                        break
                if not found_sp3:
                    continue

            # Also, check that no bond along the path is amide-like.
            peptide_bond_found = False
            for i in range(len(path) - 1):
                a1 = mol.GetAtomWithIdx(path[i])
                a2 = mol.GetAtomWithIdx(path[i+1])
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "N"}:
                    # Identify the carbon atom.
                    c_atom = a1 if a1.GetSymbol() == "C" else a2
                    double_o = False
                    has_oh = False
                    for nbr in c_atom.GetNeighbors():
                        if nbr.GetSymbol() != "O":
                            continue
                        bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                        if bond is None:
                            continue
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            double_o = True
                        elif bond.GetBondType() == Chem.BondType.SINGLE:
                            if nbr.GetTotalNumHs() >= 1 or nbr.GetFormalCharge() < 0:
                                has_oh = True
                    if double_o and not has_oh:
                        peptide_bond_found = True
                        break
            if peptide_bond_found:
                continue

            return True, f"Found free carboxylic acid and free amino group with a bond distance of {distance} (within cutoff)."

    # If no acceptable acid–amine connectivity is found, return failure.
    return False, "Functional groups found but none are close enough (or belong to the same amino acid residue) without peptide linkages."

 
# Example test cases (for development; remove or comment out in production)
if __name__ == "__main__":
    test_cases = {
        "(S)-gabaculine": "C1=C(C[C@@H](C=C1)N)C(O)=O",
        "porphyra-334": "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",
        "L-thyroxine": "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",
        "robenacoxib": "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        # False positive (should be rejected)
        "(R)-S-lactoylglutathione": "C[C@@H](O)C(=O)SC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(O)=O",
        # False negatives (expected: no matching alpha motif because free amine is not found)
        "N-Methylanthranilamide": "CNC(=O)c1ccccc1N",
        "N-(methoxyacetyl)-4-hydroxyproline": "COCC(=O)N1CC(O)CC1C(O)=O",
    }
    for name, sm in test_cases.items():
        result, reason = is_amino_acid(sm)
        print(f"{name}: {result} ({reason})")