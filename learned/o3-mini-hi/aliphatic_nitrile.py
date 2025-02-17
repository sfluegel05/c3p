"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile 
Definition: Any nitrile derived from an aliphatic compound.
Heuristic:
  1. Look for nitrile groups using SMARTS: a non‐aromatic carbon (X2) triple‐bonded to a nitrogen.
  2. For each nitrile match:
     (a) Check that the nitrile carbon is non‐aromatic.
     (b) Identify its unique substituent “R” (the neighbor that is not the nitrile nitrogen) and require that it is a carbon.
  3. From R, perform a depth‐first search (DFS) to “collect” the contiguous branch of atoms.
     Only traverse atoms that are non‐aromatic (thus allowing aliphatic unsaturation) so that we remain “in” the aliphatic part.
  4. Compute two properties of the branch:
       • carbon_fraction = (# carbon atoms)/(# heavy atoms in branch)
       • isolation: no atom in the branch should be directly attached to an atom outside the branch that is aromatic.
  5. Also, if any atom in the branch is an oxygen that is double‐bonded to a carbon (i.e. a carbonyl oxygen), then disqualify.
  6. If the branch satisfies these criteria (we require in our heuristic a carbon fraction ≥0.6),
     we consider the nitrile “attached exclusively to an aliphatic substituent branch”.
  7. If any nitrile in the molecule qualifies, return True plus an appropriate message.
     
Note: This heuristic is imperfect and may miss borderline cases.
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.

    The algorithm:
      1. Parses the SMILES and locates nitrile groups (non‐aromatic carbon triple-bonded to a nitrogen).
      2. For each nitrile group, it identifies the substituent branch (R) attached to the nitrile carbon.
      3. It performs a DFS from that R atom (only walking over non‐aromatic atoms) to collect the branch.
      4. It then computes two properties:
          - The fraction of carbons in the branch; we require at least 60% carbons.
          - None of the branch’s atoms should be directly connected to an external aromatic atom.
          - Also, if any branch atom is an oxygen in a double bond (a carbonyl), the branch is rejected.
      5. If any nitrile group passes these tests, the function returns True and a message.
         Otherwise, it returns False along with a brief explanation.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule contains a qualifying aliphatic nitrile group, False otherwise.
        str: Explanation text.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a nitrile group: a non‐aromatic carbon (with exactly 2 neighbors) triple-bonded to a nitrogen.
    nitrile_pattern = Chem.MolFromSmarts("[C;!a;X2]#[N;X1]")
    if nitrile_pattern is None:
        return False, "Error creating nitrile SMARTS pattern"

    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile group found in the molecule"

    # Helper to check if an oxygen atom has a double bond to a carbon (i.e. is in a carbonyl group)
    def is_carbonyl_oxygen(atom):
        if atom.GetSymbol() != "O":
            return False
        for bond in atom.GetBonds():
            # Check if bond is double and other atom is carbon.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetSymbol() == "C":
                    return True
        return False

    # For each nitrile group found…
    for match in matches:
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # (A) Ensure the nitrile carbon is non‐aromatic.
        if nitrile_c.GetIsAromatic():
            continue

        # (B) Identify the unique substituent R (neighbor that is not the nitrile nitrogen).
        neighbors = [nbr for nbr in nitrile_c.GetNeighbors() if nbr.GetIdx() != nitrile_n.GetIdx()]
        if not neighbors:
            continue  # unusual case: nitrile carbon with no other substituent
        R = neighbors[0]
        if R.GetSymbol() != "C":
            continue  # require that the substituent is a carbon
        
        # (C) Do a DFS from R to collect all connected atoms in the branch.
        # We only traverse atoms that are non-aromatic. (This allows sp2 unsaturated, if not aromatic.)
        branch_atom_indices = set()
        queue = [R.GetIdx()]
        while queue:
            idx = queue.pop(0)
            if idx in branch_atom_indices:
                continue
            branch_atom_indices.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Do not go back to nitrile carbon.
                if nbr.GetIdx() == nitrile_c.GetIdx():
                    continue
                # Only continue if the neighbor is non-aromatic.
                if not nbr.GetIsAromatic():
                    queue.append(nbr.GetIdx())
        if not branch_atom_indices:
            continue

        # (D) Compute the "carbon fraction" in the branch.
        total_heavy = 0
        carbon_count = 0
        for idx in branch_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() > 1:
                total_heavy += 1
                if atom.GetSymbol() == "C":
                    carbon_count += 1
        # Avoid division by zero; require branch to have at least one heavy atom.
        if total_heavy == 0:
            continue

        carbon_fraction = carbon_count / total_heavy

        # (E) Check that none of the branch atoms is an oxygen in a carbonyl group.
        branch_has_carbonyl = False
        for idx in branch_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if is_carbonyl_oxygen(atom):
                branch_has_carbonyl = True
                break
        if branch_has_carbonyl:
            continue

        # (F) Check isolation: none of the branch atoms should have a neighbor (outside the branch)
        # that is aromatic. This helps to rule out cases where the branch is directly connected to an aromatic unit.
        branch_is_isolated = True
        for idx in branch_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in branch_atom_indices:
                    if nbr.GetIsAromatic():
                        branch_is_isolated = False
                        break
            if not branch_is_isolated:
                break
        
        # (G) Now apply our heuristic:
        # We require that the branch be mostly aliphatic (at least 60% carbon) and be isolated from aromatic moieties.
        if carbon_fraction < 0.6:
            continue
        if not branch_is_isolated:
            continue

        # If we reach here, this nitrile qualifies.
        return True, "Contains a nitrile group attached to an exclusively aliphatic substituent branch"

    return False, "Nitrile group(s) found but none have a qualifying aliphatic substituent branch"


# Example usages (for testing – uncomment as needed):
# print(is_aliphatic_nitrile("N#CCC#N"))          # malononitrile -> expected True.
# print(is_aliphatic_nitrile("OC(=O)CNCC#N"))      # N-(cyanomethyl)glycine -> expected True.
# print(is_aliphatic_nitrile("CCN(CC)C(=O)CC#N"))   # N,N-diethylcyanoacetamide -> expected True.
# print(is_aliphatic_nitrile("C[C@H](N)C#N"))       # (S)-alpha-aminopropionitrile -> expected True.
# print(is_aliphatic_nitrile("O1C2=C(C(C3CCC=CC3)C(=C1N)C#N)C(=O)CC(C2)(C)C"))  # example false positive