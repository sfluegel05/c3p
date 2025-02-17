"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
Definition: An alkylamine in which one or more alkyl substituents bears an aromatic group.
Operationally, we require that at least one non‐aromatic amine nitrogen (i.e. not in any aromatic ring)
has at least one carbon-only substituent that (via a chain of single C–C bonds) leads to an aromatic carbon.
This chain “branch” is assumed to be an alkyl chain. To avoid spurious results we only traverse C atoms,
ignoring any branch that goes off into carbonyls or heteroatoms.
"""

from rdkit import Chem
from collections import deque

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    The approach is:
      1. Parse the SMILES.
      2. Identify non‐aromatic amine nitrogen atoms (atomic number 7 not flagged as aromatic).
      3. For each such nitrogen, check for any substituent branch that is attached via a SINGLE bond 
         to a carbon atom. If that immediate neighbor is aromatic then report immediately.
      4. Otherwise, follow that branch via a breadth‐first search (BFS) restricted to carbon atoms and only
         following single bonds. If an aromatic carbon is encountered within a maximum total bond distance of 3
         (i.e. including the bond from the N) then report success.
      5. Otherwise, if no branch of any non‐aromatic amine nitrogen meets these criteria, return False.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      (bool, str): A tuple of a boolean indicating whether the molecule qualifies as an aralkylamine,
                   plus a string explanation.
    
    If the SMILES cannot be parsed, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    max_depth = 3  # Maximum bond distance (including the immediate bond from nitrogen)

    # Iterate over all atoms. Consider only nitrogen atoms (atomic number 7)
    # that are not flagged as aromatic.
    for n_atom in mol.GetAtoms():
        if n_atom.GetAtomicNum() != 7 or n_atom.GetIsAromatic():
            continue  # Skip if not a non-aromatic nitrogen.
        
        n_idx = n_atom.GetIdx()
        # Examine each bond from the nitrogen.
        for bond in n_atom.GetBonds():
            # Only consider single bonds (alkyl-type bond)
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            neighbor = bond.GetOtherAtom(n_atom)
            # We only follow branches that start with a carbon.
            if neighbor.GetAtomicNum() != 6:
                continue

            # For the branch attached at nitrogen, get the immediate neighbor index.
            start_idx = neighbor.GetIdx()
            # If the immediate neighbor is aromatic then it is already directly an aromatic substituent.
            if neighbor.GetIsAromatic():
                reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) directly attached to an "
                          f"aromatic carbon (atom {start_idx}) via a single bond (distance 1). "
                          "Molecule classified as aralkylamine.")
                return True, reason

            # Otherwise, perform a BFS along the branch; we restrict allowed atoms to carbon atoms 
            # connected by single bonds.
            # The BFS queue stores tuples (current_atom_idx, distance) where distance includes the bond from N.
            visited = set([n_idx, start_idx])
            queue = deque([(start_idx, 1)])  # starts at distance 1 from the nitrogen

            while queue:
                current_idx, depth = queue.popleft()
                # If we reached maximum allowed distance, do not continue beyond.
                if depth >= max_depth:
                    # But still check the current node in case it is aromatic.
                    current_atom = mol.GetAtomWithIdx(current_idx)
                    if current_atom.GetAtomicNum() == 6 and current_atom.GetIsAromatic():
                        reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) with an aromatic substituent "
                                  f"(atom {current_idx}) at a bond distance of {depth}. Molecule classified as aralkylamine.")
                        return True, reason
                    continue

                # Expand the search: only allowed atoms are carbons connected by single bonds.
                current_atom = mol.GetAtomWithIdx(current_idx)
                for nbr in current_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in visited:
                        continue
                    bond_ = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
                    if bond_.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue
                    # Only traverse if the atom is carbon.
                    if nbr.GetAtomicNum() != 6:
                        continue
                    new_depth = depth + 1
                    # If nbr is aromatic and new_depth is within allowed limit then we have an eligible branch.
                    if nbr.GetIsAromatic():
                        reason = (f"Found a non‐aromatic amine nitrogen (atom {n_idx}) with an aromatic substituent "
                                  f"(atom {nbr_idx}) at a bond distance of {new_depth}. Molecule classified as aralkylamine.")
                        return True, reason
                    # Otherwise add the neighbor to the search if we are not past the limit.
                    if new_depth < max_depth:
                        visited.add(nbr_idx)
                        queue.append((nbr_idx, new_depth))
            # End BFS for this branch.
    # If no branch resulted in an aromatic carbon found through a pure alkyl (C-only) chain, then fail.
    return False, ("No aralkylamine substructure found: "
                   "no non‐aromatic amine nitrogen is connected via an alkyl chain (single-bonded carbon branch) "
                   "to an aromatic carbon within 3 bonds.")

# Optional test block; remove or comment out when deploying as a module.
if __name__ == "__main__":
    # Some example test cases.
    test_cases = [
        ("NCc1ccccc1", "benzylamine"),  # True, branch: N-CH2 (aliphatic) then aromatic ring
        ("NCCc1ccccc1", "phenethylamine"),  # True, branch length 2
        ("c1ccc(N)cc1", "aniline (should not classify, N aromatic)"),  # False: nitrogen is aromatic
        ("OCCNC1=CC=CC=C1", "2-Anilinoethanol"),  # True; branch from N to CH2 to aromatic ring
    ]
    
    for smi, name in test_cases:
        res, explanation = is_aralkylamine(smi)
        print(f"SMILES: {smi}\nMolecule: {name}\nClassification: {res}\nReason: {explanation}\n")