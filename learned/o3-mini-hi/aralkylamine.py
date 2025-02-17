"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine 
Definition: An alkylamine in which the alkyl group is substituted by an aromatic group.
That is, a non‐aromatic (aliphatic) N atom attached (through one or more sp3 carbons)
to an aromatic ring (but not directly attached to an aromatic atom).
"""

from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Approach:
      1. Parse the SMILES into an RDKit molecule.
      2. Find each nitrogen atom that is not part of an aromatic ring 
         (i.e. an aliphatic amine).
      3. For each such nitrogen, examine its neighbors. For any neighbor that is a carbon 
         and is sp3-hybridized (i.e. likely part of an alkyl chain), perform a breadth-first 
         search along the chain. If the search (starting at depth>=2, to avoid directly attached aromatic atoms)
         finds an aromatic atom, then we classify the molecule as an aralkylamine.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as an aralkylamine, False otherwise.
        str: An explanation for the classification.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper: for a starting atom, do a breadth-first search (BFS) to see if we can reach an aromatic atom.
    # We require that the aromatic atom is not directly attached to the nitrogen (i.e. depth>=2).
    def bfs_for_aromatic(start_atom, origin_atom):
        from collections import deque
        # Each item: (atom, depth, came_from)
        queue = deque()
        queue.append((start_atom, 1, origin_atom))
        visited = {start_atom.GetIdx(), origin_atom.GetIdx()}
        while queue:
            current_atom, depth, parent = queue.popleft()
            # If we are at depth>=2 and current atom is aromatic then we have a hit.
            if depth >= 2 and current_atom.GetIsAromatic():
                return True
            # Otherwise, proceed to neighbors.
            for nbr in current_atom.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx in visited:
                    continue
                # Only continue along non-aromatic carbons (or other atoms) to stay in the alkyl chain.
                # (We do not force all atoms along the chain to be carbon but one could refine this.)
                # Here we simply add the neighbor if it is not the origin (amine N) to avoid going back.
                visited.add(nidx)
                queue.append((nbr, depth+1, current_atom))
        return False

    # Look for at least one non-aromatic nitrogen that has a suitable substituent chain.
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen.
        if atom.GetAtomicNum() == 7:
            # Check if nitrogen is non-aromatic: this is our aliphatic amine.
            if atom.GetIsAromatic():
                continue  # Skip if the amine is directly aromatic.
            # Examine neighbors of the nitrogen.
            for nbr in atom.GetNeighbors():
                # We consider the substituent if it is a carbon and sp3-hybridized.
                if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    # Use BFS from this carbon. We pass the nitrogen as the origin
                    # so we can require a minimum chain length of 2 bonds from the nitrogen.
                    if bfs_for_aromatic(nbr, atom):
                        return True, ("Found an aliphatic amine (N not in aromatic ring) with an alkyl chain substituent " 
                                      "that leads to an aromatic group. Molecule classified as aralkylamine.")
    return False, "No aralkylamine substructure found: no aliphatic amine found with an alkyl chain reaching an aromatic group."

# Example test cases (can be removed or commented out when using as a module)
if __name__ == "__main__":
    # benzylamine should be aralkylamine: NH-CH2-C6H5
    smiles_examples = [
        "NCc1ccccc1",  # benzylamine
        "NCCc1ccccc1",  # phenethylamine
        "c1ccc(N)cc1"   # aniline (should not count, as amine directly attached to aromatic ring)
    ]
    for sm in smiles_examples:
        result, reason = is_aralkylamine(sm)
        print(f"SMILES: {sm}\nClassification: {result}\nReason: {reason}\n")