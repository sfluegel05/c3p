"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: Sulfolipid
A sulfolipid is defined as a compound containing a sulfonic acid residue 
(joined either directly (C–S) or via one bridging oxygen (C–O–S)) 
that is connected to a lipid moiety. The lipid is defined here as having 
a long, largely unbranched, aliphatic chain (of at least 12 sp³ carbons) 
and a degree of molecular complexity (e.g. presence of at least one ring).
Additionally, the molecular weight must exceed 300 Da.
This heuristic aims to minimize false positives such as simple alkyl sulfonates.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    
    The criteria are:
      - Molecular weight > 300 Da.
      - Contains a sulfonate/sulfonic acid group (S(=O)(=O)[O-] or S(=O)(=O)O).
      - Contains a long, essentially unbranched, aliphatic chain of at least 12 sp3 carbons.
      - The sulfonate must be “linked” (directly or via a single oxygen) 
        to this long chain (i.e. within 3 bonds).
      - The overall molecular complexity (e.g. presence of rings) must suggest
        that the molecule is not a mere simple alkyl sulfonate.
    
    Args:
       smiles (str): The SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule is classified as a sulfolipid, False otherwise.
       str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a sulfolipid (mol wt = {mol_wt:.1f})"
    
    # Define sulfonate SMARTS patterns: charged or neutral sulfonic acid
    sulfo_smarts_list = [ "S(=O)(=O)[O-]", "S(=O)(=O)O" ]
    sulfonate_idx_set = set()
    for smarts in sulfo_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(patt, useChirality=True)
        for match in matches:
            # assume the sulfur is the first atom in our SMARTS.
            sulfonate_idx_set.add(match[0])
    if not sulfonate_idx_set:
        return False, "No sulfonic acid group detected"
    
    # Identify molecular complexity: count rings.
    rings = mol.GetRingInfo().NumRings()
    
    # Helper function: compute longest linear (unbranched) alkyl chain.
    # We consider eligible atoms as aliphatic (sp3, non‐aromatic, non‐ring) carbons.
    def get_longest_aliphatic_chain(mol):
        eligible = []
        for atom in mol.GetAtoms():
            # Identify carbons that are sp3, not aromatic, and not in any ring.
            if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                if not atom.GetIsAromatic() and not atom.IsInRing():
                    eligible.append(atom.GetIdx())
        if not eligible:
            return 0, []
        # Build a simple graph over eligible atoms
        # We use a brute force approach: for each eligible atom, do a BFS to find the longest path
        best_length = 0
        best_path = []
        for start in eligible:
            # visited: atom idx -> path (list of idxs) from start
            visited = {start: [start]}
            queue = [start]
            while queue:
                current = queue.pop(0)
                for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in eligible and nbr_idx not in visited:
                        path = visited[current] + [nbr_idx]
                        visited[nbr_idx] = path
                        queue.append(nbr_idx)
                        if len(path) > best_length:
                            best_length = len(path)
                            best_path = path
        return best_length, best_path

    longest_chain_length, longest_chain_atoms = get_longest_aliphatic_chain(mol)
    if longest_chain_length < 12:
        return False, f"No sufficiently long aliphatic chain found (maximum chain length = {longest_chain_length})"
    
    # Now, assess connectivity between a sulfonate group and the lipid chain.
    # The sulfonate S must be connected via 1 or 2 bonds (i.e. path length <= 3)
    found_connection = False
    connection_reason = ""
    for s_idx in sulfonate_idx_set:
        for atom_idx in longest_chain_atoms:
            path = rdmolops.GetShortestPath(mol, s_idx, atom_idx)
            # len(path)-1 gives number of bonds in the shortest route.
            if (len(path) - 1) <= 3:
                found_connection = True
                connection_reason = (f"Found sulfonic acid group (S at atom idx {s_idx}) connected "
                                     f"to a long aliphatic chain (length = {longest_chain_length}) via a path "
                                     f"of {len(path)-1} bonds.")
                break
        if found_connection:
            break
    if not found_connection:
        return False, "Sulfonic acid group not connected to a long aliphatic chain within 3 bonds"
    
    # Reject overly simple molecules that appear to be nothing but an alkyl sulfonate.
    # Here we consider a lack of rings (and other complexity) and an overall heavy atom count
    # that is very similar to that of the long chain plus the sulfonate group.
    num_heavy = mol.GetNumHeavyAtoms()
    # Assume roughly 4 atoms in sulfonate (+ connector); if the heavy atoms are only the long chain and sulfonate,
    # then we consider it too simple.
    if rings == 0 and num_heavy <= (longest_chain_length + 4):
        return False, "Molecule appears to be a simple alkyl sulfonate, not a sulfolipid"
    
    return True, connection_reason

# Example usage (for testing):
if __name__ == '__main__':
    test_smiles = [
        # psychosine sulfate: should be a sulfolipid.
        "CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](COS(O)(=O)=O)[C@H](O)[C@H](O)[C@H]1O",
        # (3'-SulfO)galbeta-cer(D18:1/20:0): should be a sulfolipid.
        "S(OC1[C@@H](O)[C@H](O[C@@H](OC[C@H](NC(=O)CCCCCCCCCCCCCCCCCCC)[C@H](O)/C=C/CCCCCCCCCCCCC)C1O)CO)(O)(=O)=O",
        # hexadecane-1-sulfonate (a false positive previously): should not be a sulfolipid.
        "CCCCCCCCCCCCCCCS(O)(=O)=O",
    ]
    for smi in test_smiles:
        is_sulfo, reason = is_sulfolipid(smi)
        print("SMILES:", smi)
        print("Classified as sulfolipid?", is_sulfo)
        print("Reason:", reason)
        print("-----")