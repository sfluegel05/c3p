"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies chemical entities of the class amino acid.
Definition used: An amino acid is a molecule containing at least one free (non‑acylated, primary, sp³) amino group and at least one carboxylic acid group 
(in either –C(=O)[OH] or –C(=O)[O-] form). In addition, an acid–amine pair must appear “close together” (ideally sharing an “alpha‐carbon”,
i.e. the acid carbon is directly bonded to an sp3 carbon that in turn is bonded to a free amino nitrogen). 
Some borderline cases may remain ambiguous.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is (likely) an amino acid based on its SMILES string.

    The algorithm is as follows:
      1. Parse the molecule.
      2. Identify carboxylic acid groups (both neutral and anionic) where the acid carbon 
         (the carbonyl carbon) has degree = 3.
      3. Identify free (non‑acylated) amino groups. A free amino group in our criteria must:
            a. Be a nitrogen having formal charge of zero.
            b. Not be aromatic.
            c. Have at least one attached hydrogen (suggesting it is primary).
            d. Not be “acylated” (i.e. bonded via a carbon that itself bears a double‐bonded oxygen without an –OH).
      4. For each acid group, check if any neighboring carbon (candidate α‑carbon) is sp³ and is directly
         connected to one of the free amino groups. If such an “α‐motif” is found (bond distance = 2), return True.
      5. If not, then for every acid–amine pair in the same connected fragment, compute their shortest
         bond path. Accept the pair if the distance is ≤ 4, at least one intermediate atom (if any) is sp³, and
         none of the bonds along the path appear “amide‐like.” 
      6. If no acceptable close connectivity is found, return False.
     
    Args:
       smiles (str): SMILES string of the molecule.
        
    Returns:
       bool: True if the molecule is classified (likely) as an amino acid, False otherwise.
       str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ---------------- Step 1: Identify acid groups ----------------
    # We look for neutral acids and anionic acids.
    acid_neutral_smarts = "[CX3](=O)[OX2H]"
    acid_anion_smarts   = "[CX3](=O)[O-]"
    acid_neutral = Chem.MolFromSmarts(acid_neutral_smarts)
    acid_anion   = Chem.MolFromSmarts(acid_anion_smarts)
    acid_indices = set()  # will contain indices of acid carbons
    for match in mol.GetSubstructMatches(acid_neutral):
        acid_atom = mol.GetAtomWithIdx(match[0])
        if acid_atom.GetDegree() == 3:
            acid_indices.add(match[0])
    for match in mol.GetSubstructMatches(acid_anion):
        acid_atom = mol.GetAtomWithIdx(match[0])
        if acid_atom.GetDegree() == 3:
            acid_indices.add(match[0])
    if not acid_indices:
        return False, "No carboxylic acid group found"

    # ---------------- Step 2: Identify free (non‑acylated) amino groups ----------------
    free_amino_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        if atom.GetFormalCharge() != 0:
            continue
        # Reject if the nitrogen is aromatic (we expect a free aliphatic amine)
        if atom.GetIsAromatic():
            continue
        # Prefer primary amines: require at least one hydrogen.
        if atom.GetTotalNumHs() < 1:
            continue

        # Check if this N is acylated: for any neighboring carbon, if that carbon is bound to a double bonded oxygen
        # (and not to an –OH), then we consider the nitrogen acylated.
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
                        if nn.GetTotalNumHs() >= 1 or nn.GetFormalCharge() < 0:
                            has_hydroxyl = True
            if has_carbonyl and not has_hydroxyl:
                is_acylated = True
                break
        if not is_acylated:
            free_amino_indices.append(atom.GetIdx())
    if not free_amino_indices:
        return False, "No free (non‑acylated, primary) amino group found"

    # ---------------- Step 3: Look for direct α‐motif (ideal amino acid core) ----------------
    for acid_idx in acid_indices:
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            # Candidate α‑carbon should be sp³ (not aromatic) 
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Now check if this candidate has a free amino group as a neighbor (other than the acid carbon)
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == acid_idx:
                    continue
                if nbr2.GetSymbol() == "N" and nbr2.GetIdx() in free_amino_indices:
                    return True, ("Found acid carbon bonded to an sp³ alpha‑carbon which is directly bonded to a free amino group "
                                  "(bond distance = 2).")
    # ---------------- Step 4: Fallback – search via shortest path between acid and free amine ----------------
    # Only consider acid and amine that belong to the same connected fragment.
    frags = list(rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False))
    for acid_idx in acid_indices:
        for amine_idx in free_amino_indices:
            # Check they are in the same fragment.
            acid_frag = None
            amine_frag = None
            for frag in frags:
                if acid_idx in frag:
                    acid_frag = frag
                if amine_idx in frag:
                    amine_frag = frag
            if acid_frag is None or amine_frag is None or acid_frag != amine_frag:
                continue

            path = rdmolops.GetShortestPath(mol, acid_idx, amine_idx)
            if not path or len(path) < 2:
                continue
            distance = len(path) - 1  # number of bonds
            if distance > 4:
                continue  # too far apart to be within a single amino acid-like residue

            # For fallback paths longer than 2 bonds, require that at least one intermediate atom is sp³
            if distance > 2:
                has_sp3 = any(mol.GetAtomWithIdx(idx).GetHybridization() == Chem.rdchem.HybridizationType.SP3 
                              for idx in path[1:-1])
                if not has_sp3:
                    continue

            # Also make sure no bond in the path is amide-like (i.e. a C–N bond where the C is in a carbonyl lacking –OH)
            peptide_bond_found = False
            for i in range(len(path) - 1):
                a1 = mol.GetAtomWithIdx(path[i])
                a2 = mol.GetAtomWithIdx(path[i+1])
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "N"}:
                    # Identify the carbon atom in the C–N connection.
                    c_atom = a1 if a1.GetSymbol() == "C" else a2
                    double_o = False
                    has_oh = False
                    for nbr in c_atom.GetNeighbors():
                        if nbr.GetSymbol() != "O":
                            continue
                        bnd = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                        if bnd is None:
                            continue
                        if bnd.GetBondType() == Chem.BondType.DOUBLE:
                            double_o = True
                        elif bnd.GetBondType() == Chem.BondType.SINGLE:
                            if nbr.GetTotalNumHs() >= 1 or nbr.GetFormalCharge() < 0:
                                has_oh = True
                    if double_o and not has_oh:
                        peptide_bond_found = True
                        break
            if peptide_bond_found:
                continue

            return True, (f"Found free carboxylic acid and free amino group with a bond distance of {distance} "
                          "within cutoff (fallback connectivity).")
    return False, ("Functional groups found but none appear sufficiently connected (no clear acid–amine motif "
                   "within 4 bonds that avoids amide-like connectivity).")

# ---------------- Optional test cases (remove or comment out in production) ----------------
if __name__ == "__main__":
    # A selection of examples from the outcome summary.
    test_cases = {
        "(S)-gabaculine": "C1=C(C[C@@H](C=C1)N)C(O)=O",
        "porphyra-334": "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",
        "L-thyroxine": "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",
        "robenacoxib": "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        "(R)-S-lactoylglutathione": "C[C@@H](O)C(=O)SC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(O)=O",
        "N-Methylanthranilamide": "CNC(=O)c1ccccc1N",
        "N-(methoxyacetyl)-4-hydroxyproline": "COCC(=O)N1CC(O)CC1C(O)=O",
    }
    for name, sm in test_cases.items():
        result, reason = is_amino_acid(sm)
        print(f"{name}: {result} ({reason})")