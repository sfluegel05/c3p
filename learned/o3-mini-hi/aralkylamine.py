"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
Definition: An alkylamine in which one or more alkyl substituents bears an aromatic group.
Operationally, we require that at least one non‐aromatic amine nitrogen 
(i.e. a nitrogen not flagged as aromatic) is connected — directly or via an alkyl chain (of 2–3 bonds)
that is composed solely of saturated (sp3) carbon atoms—to an aromatic carbon.
Direct attachment of an aromatic carbon (i.e. an arylamine) is not acceptable.
For cyclic amines the search also examines substituents coming off ring atoms that are attached 
to the amine.
"""

from rdkit import Chem
from rdkit.Chem import rdchem
from collections import deque

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Approach:
      1. Parse the SMILES.
      2. For each nitrogen (atomic num 7) that is non-aromatic, examine its substituents.
         (Note: direct attachment of an aromatic carbon is not allowed, because that defines arylamine.)
      3. For substituents that start with a carbon that is sp3 and non–aromatic, do a breadth-first
         search limited to a maximum total bond distance of 3 from the nitrogen.
         Only single bonds are traversed, and only carbon atoms (atomic num 6) with sp3 hybridization
         and not aromatic are allowed along the branch.
         If an aromatic carbon is encountered (as a terminal hit) with a bond distance of 2 or 3, 
         the aralkylamine criterion is satisfied.
      4. For cyclic amines the search is extended: if the amine is in a ring then for each ring neighbor,
         we check its exocyclic substituents (neighbors not in the ring) as potential branch starts.
      5. If any branch meets the criteria, returns True and a message; otherwise, returns False.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple with boolean classification and a message explaining the decision.
                   If SMILES cannot be parsed, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    max_depth = 3  # maximum allowed bond distance from the nitrogen
    # Helper: BFS search along an alkyl branch.
    def bfs_branch(start_atom, start_depth, visited):
        # Use a deque for BFS; each entry is a tuple (atom, current depth)
        queue = deque([(start_atom, start_depth)])
        while queue:
            current_atom, depth = queue.popleft()
            # Check if we've reached an aromatic carbon (terminal hit) at valid branch length:
            if depth >= 2 and current_atom.GetIsAromatic():
                return current_atom.GetIdx(), depth
            # Stop if we reached our maximum allowed depth.
            if depth >= max_depth:
                continue
            # Expand: traverse only single bonds to carbon atoms that are sp3 and non-aromatic.
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Skip if already visited in this branch.
                if nbr_idx in visited:
                    continue
                bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), nbr_idx)
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Two cases: if neighbor is aromatic, then if reached at valid depth we can immediately accept.
                if nbr.GetIsAromatic():
                    if (depth + 1) >= 2 and (depth + 1) <= max_depth:
                        return nbr.GetIdx(), depth + 1
                    continue
                # Must be a carbon (atomic number 6) that is sp3 and non-aromatic.
                if nbr.GetAtomicNum() != 6:
                    continue
                if nbr.GetHybridization() != rdchem.HybridizationType.SP3:
                    continue
                visited.add(nbr_idx)
                queue.append((nbr, depth + 1))
        return None

    # Now iterate over all atoms to find a suitable nitrogen.
    for n_atom in mol.GetAtoms():
        if n_atom.GetAtomicNum() != 7:
            continue
        if n_atom.GetIsAromatic():
            continue  # skip aromatic nitrogens
        n_idx = n_atom.GetIdx()
        # First, consider all substituents directly attached to the N that are not aromatic.
        for nbr in n_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(n_idx, nbr.GetIdx())
            if bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            # We want the branch to start with a carbon.
            if nbr.GetAtomicNum() != 6:
                continue
            # If the carbon is directly aromatic, that would be an arylamine – skip.
            if nbr.GetIsAromatic():
                continue
            # Require that the first carbon is sp3.
            if nbr.GetHybridization() != rdchem.HybridizationType.SP3:
                continue
            # For nitrogens that are not in a ring or for substituents that are exocyclic,
            # do a BFS starting with depth = 1.
            if not n_atom.IsInRing() or (not nbr.IsInRing()):
                # Create a new visited set with the starting atoms.
                visited = {n_idx, nbr.GetIdx()}
                res = bfs_branch(nbr, 1, visited)
                if res:
                    target_idx, used_depth = res
                    reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) with an aromatic substituent "
                              f"(atom {target_idx}) at a bond distance of {used_depth}. Molecule classified as aralkylamine.")
                    return True, reason
            # Now, if the nitrogen is in a ring, it is common that its immediate neighbors are ring atoms.
            # In that case, check for exocyclic branches from the ring neighbor.
            if n_atom.IsInRing() and nbr.IsInRing():
                # For each neighbor of this ring atom (nbr) that is not in the ring, consider it a branch start.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == n_idx:
                        continue
                    # We want an exocyclic branch: nbr2 should not be in a ring.
                    if nbr2.IsInRing():
                        continue
                    if nbr2.GetAtomicNum() != 6:
                        continue
                    if nbr2.GetIsAromatic():
                        continue  # if directly aromatic, that would be a direct aryl attachment – not allowed.
                    if nbr2.GetHybridization() != rdchem.HybridizationType.SP3:
                        continue
                    # For these branches, the bond path is: N -(bond already traversed)-> nbr (ring atom) -(bond)-> nbr2.
                    # This counts as depth = 2.
                    visited = {n_idx, nbr.GetIdx(), nbr2.GetIdx()}
                    res = bfs_branch(nbr2, 2, visited)
                    if res:
                        target_idx, used_depth = res
                        reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) via ring‐branching "
                                  f"with an aromatic substituent (atom {target_idx}) at a bond distance of {used_depth} (path including exocyclic branch). "
                                  f"Molecule classified as aralkylamine.")
                        return True, reason
    # If nothing is found:
    return (False, "No aralkylamine substructure found: "
            "no non‐aromatic amine nitrogen is connected via a suitable alkyl chain "
            "to an aromatic carbon within allowed bond distance.")

# Optional test block; remove or comment this block out when deploying as a module.
if __name__ == "__main__":
    # Some test cases including those expected to be true positives and a few negatives.
    test_cases = [
        ("NCc1ccccc1", "benzylamine (should classify True)"),
        ("NCCc1ccccc1", "phenethylamine (should classify True)"),
        ("c1ccc(N)cc1", "aniline (should classify False: N is aromatic)"),
        ("OCCNC1=CC=CC=C1", "2-Anilinoethanol (should classify True)"),
        # Trihexyphenidyl hydrochloride is a known false negative in the earlier version.
        ("Cl.C1=CC=CC(=C1)C(C2CCCCC2)(CCN3CCCCC3)O", "Trihexyphenidyl hydrochloride (should classify True)"),
    ]
    
    for smi, name in test_cases:
        res, explanation = is_aralkylamine(smi)
        print(f"SMILES: {smi}\nMolecule: {name}\nClassification: {res}\nReason: {explanation}\n")