"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies chemical entities of the class amino acid.
Definition used: An amino acid is a molecule containing at least one free (non‐acylated) amino group and at least one carboxylic acid
group (neutral –C(=O)[OH] or anionic –C(=O)[O-]). In addition, at least one acid–amine pair must appear “close together”
(i.e. with a shortest bond path of 4 bonds or fewer) and, ideally, these groups share an “alpha‐carbon” (or are directly connected).
This is a heuristic meant to capture a backbone-like motif (or its close analogue) rather than merely finding the two isolated groups.
Note: Some borderline cases may still be ambiguous.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is (likely) an amino acid based on its SMILES string.

    The algorithm requires:
       1. Valid parsing.
       2. At least one carboxylic acid group (in neutral –C(=O)[OH] or anionic –C(=O)[O-] form) whose carbon has
          the expected connectivity (degree exactly 3).
       3. At least one free (non‑acylated) amino group (where the nitrogen has formal charge 0 and at most 2 heavy‐atom neighbors).
       4. The existence of at least one acid–amine pair with a shortest bond path of 4 bonds or fewer.
          In addition, if the path is longer than 2 bonds we check that one of the intermediate atoms (a candidate “alpha‐carbon”)
          is sp³ – as expected in a real amino acid residue – and that none of the bonds along the path resemble peptide (amide) bonds.
    
    Args:
       smiles (str): SMILES string of the molecule.
        
    Returns:
       bool: True if the molecule is classified as an amino acid, False otherwise.
       str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # -------------------------
    # 1. Identify carboxyl groups.
    # We use two SMARTS patterns: one for neutral acids and one for anions.
    acid_neutral_smarts = "[CX3](=O)[OX2H]"
    acid_anion_smarts   = "[CX3](=O)[O-]"
    acid_neutral = Chem.MolFromSmarts(acid_neutral_smarts)
    acid_anion   = Chem.MolFromSmarts(acid_anion_smarts)
    acid_indices = set()  # store the carbon index for each acid group found

    for match in mol.GetSubstructMatches(acid_neutral):
        # Expect match to have (C, O_H) so that the "acid carbon" is match[0].
        acid_atom = mol.GetAtomWithIdx(match[0])
        # Check that this acid carbon has degree 3 (i.e. is bonded to the carbonyl O and the -OH, plus one C)
        if acid_atom.GetDegree() == 3:
            acid_indices.add(match[0])
    for match in mol.GetSubstructMatches(acid_anion):
        acid_atom = mol.GetAtomWithIdx(match[0])
        if acid_atom.GetDegree() == 3:
            acid_indices.add(match[0])

    if not acid_indices:
        return False, "No carboxylic acid group found"

    # -------------------------
    # 2. Identify free (non-acylated) amino groups.
    free_amino_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        # We require the nitrogen to have zero formal charge.
        if atom.GetFormalCharge() != 0:
            continue
        # Count heavy-atom neighbors (i.e. non-hydrogen)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) > 2:
            # Too heavily substituted to be a simple free amine.
            continue

        n_idx = atom.GetIdx()
        is_acylated = False
        # Check if this N is directly bonded to a carbonyl that is not part of a proper acid group.
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            # For each bond from this N–C neighbor, check if the C is bonded to a double‐bonded O but not also an –OH.
            double_o = False
            has_oh = False
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == n_idx:
                    continue
                if nbr2.GetSymbol() == "O":
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                    if bond is None:
                        continue
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_o = True
                    elif bond.GetBondType() == Chem.BondType.SINGLE:
                        if nbr2.GetTotalNumHs() >= 1 or nbr2.GetFormalCharge() < 0:
                            has_oh = True
            # If it’s an amide-style carbonyl (i.e. double-bonded O but no -OH) and not one of our acid groups, mark as acylated.
            if double_o and not has_oh:
                is_acylated = True
                break
        if not is_acylated:
            free_amino_indices.append(n_idx)

    if not free_amino_indices:
        return False, "No free (non-acylated) amino group found"

    # -------------------------
    # 3. Look for an acid–amine pair that are “close enough”
    # We restrict our attention to pairs in the same connected fragment.
    frags = list(rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False))
    for acid_idx in acid_indices:
        for amine_idx in free_amino_indices:
            # Determine if acid and amine are in same fragment.
            acid_frag = None
            amine_frag = None
            for frag in frags:
                if acid_idx in frag:
                    acid_frag = frag
                if amine_idx in frag:
                    amine_frag = frag
            if acid_frag is None or amine_frag is None or acid_frag != amine_frag:
                continue

            # Compute shortest path between acid carbon and amino nitrogen.
            path = rdmolops.GetShortestPath(mol, acid_idx, amine_idx)
            if not path or len(path) < 2:
                continue
            distance = len(path) - 1  # number of bonds
            if distance > 4:
                continue  # too long to be part of the same residue

            # Extra heuristic: if distance >= 3, try to see if one intermediate atom appears to be an sp3
            # "connecting" (alpha) carbon.
            if distance >= 3:
                # For path lengths 3 or 4, check the atoms in the middle (excluding acid C and amine N).
                candidate_alpha = None
                for idx in path[1:-1]:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        candidate_alpha = atom
                        break
                if candidate_alpha is None:
                    # If none of the intermediate atoms is sp3, then the groups may be too remote (or overlapping aromatic systems),
                    # so we skip this pair.
                    continue

            # -------------------------
            # 4. Inspect bonds along path to ensure no intervening peptide (amide) bond.
            peptide_bond_found = False
            for i in range(len(path)-1):
                a1 = mol.GetAtomWithIdx(path[i])
                a2 = mol.GetAtomWithIdx(path[i+1])
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "N"}:
                    # Identify the C atom.
                    c_atom = a1 if a1.GetSymbol() == "C" else a2
                    double_o = False
                    has_oh = False
                    for nbr in c_atom.GetNeighbors():
                        # Skip if the neighbor is not oxygen.
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
                        # This bond is amide-like.
                        peptide_bond_found = True
                        break
            if peptide_bond_found:
                continue

            # If we have found an acid–amine pair with a short bond path and minimal acylation, we classify as amino acid.
            return True, f"Found free carboxylic acid and free amino group with a bond distance of {distance} (within cutoff)."

    # If no acceptable acid–amine pair is found then reject.
    return False, "Functional groups found but none are close enough (or belong to the same amino acid residue) without peptide linkages."

# Example test cases (for development; remove or comment out in production)
if __name__ == "__main__":
    test_cases = {
        "(S)-gabaculine": "C1=C(C[C@@H](C=C1)N)C(O)=O",
        "porphyra-334": "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",
        "L-thyroxine": "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",
        "robenacoxib": "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        # Some compounds that were false positives previously:
        "O-linoelaidylcarnitine": "CCCCC\\C=C\\C\\C=C\\CCCCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C",
        "(R)-S-lactoylglutathione": "C[C@@H](O)C(=O)SC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(O)=O",
        # Some false negatives previously:
        "N-Methylanthranilamide": "CNC(=O)c1ccccc1N",
        "N-(methoxyacetyl)-4-hydroxyproline": "COCC(=O)N1CC(O)CC1C(O)=O",
    }
    for name, sm in test_cases.items():
        result, reason = is_amino_acid(sm)
        print(f"{name}: {result} ({reason})")