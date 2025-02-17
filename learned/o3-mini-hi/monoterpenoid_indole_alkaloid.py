"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid

A monoterpenoid indole alkaloid is biosynthesised from L-tryptophan and a diisoprenoid (usually secologanin)
building block. The heuristic approach:
  1. Identify an indole core (free, N-substituted, or via fused ring system).
  2. Look for a non‐aromatic C=C bond near the indole core.
     Then “grow” the attached fragment (excluding the indole) and require that this fragment contains
     at least a threshold number of carbon atoms (proxy for a terpenoid moiety).
  3. Basic size check (total carbon count and molecular weight).
  
Note: This heuristic will not capture every nuance of biosynthesis.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid using heuristic substructure patterns.
    
    The classification is based on:
      - The presence of an indole core (free, substituted, or detected via ring fusion).
      - The presence of a non‐aromatic vinyl (C=C) bond close (within 4 bonds) to the indole core.
        Moreover, the vinyl bond should be part of a fragment (excluding the indole core)
        that contains at least 5 carbon atoms (proxy for a terpenoid unit).
      - Basic size requirements in terms of carbon count and molecular weight.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoterpenoid indole alkaloid, False otherwise.
        str: A reason describing the classification outcome.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- 1. Identify the indole core -----
    # Try SMARTS for free indole and N-substituted indole.
    indole_free = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c2")
    indole_sub  = Chem.MolFromSmarts("c1ccc2c(c1)[n]c2")
    
    indole_atoms = None
    if mol.HasSubstructMatch(indole_free):
        indole_atoms = set(mol.GetSubstructMatch(indole_free))
    elif mol.HasSubstructMatch(indole_sub):
        indole_atoms = set(mol.GetSubstructMatch(indole_sub))
    else:
        # If simple SMARTS did not match, try to find a fused ring system.
        rings = mol.GetRingInfo().AtomRings()
        for ring1 in rings:
            if len(ring1) == 5:
                # Check if the 5-membered ring has at least one nitrogen
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring1):
                    for ring2 in rings:
                        if len(ring2) == 6:
                            # Check that every atom in the 6-membered ring is aromatic.
                            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring2):
                                # Must share at least 2 atoms for a fused system.
                                if len(set(ring1).intersection(ring2)) >= 2:
                                    indole_atoms = set(ring1).union(ring2)
                                    break
                    if indole_atoms is not None:
                        break
    if indole_atoms is None:
        return False, "No indole core (free, substituted, or fused) found"
    
    # ----- 2. Look for an attached terpene-like fragment -----
    # We search for non-aromatic C=C bonds.
    terpene_found = False
    dm = Chem.GetDistanceMatrix(mol)  # distance matrix for path lengths
    
    # Helper function: given a starting atom idx, perform a BFS over atoms that are NOT in the indole core.
    def bfs_fragment(start_idx, excluded_set):
        visited = set()
        queue = [start_idx]
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                nid = neighbor.GetIdx()
                # Only traverse atoms not in the indole core.
                if nid not in excluded_set and nid not in visited:
                    queue.append(nid)
        return visited

    # Iterate over bonds to find non-aromatic double bonds between carbons.
    candidate_bonds = []
    for bond in mol.GetBonds():
        # Check if bond is a double bond (using GetBondTypeAsDouble()==2) and is non-aromatic.
        if bond.GetBondTypeAsDouble() == 2 and not bond.GetIsAromatic():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)
            # Only consider if both atoms are carbons.
            if atom_i.GetAtomicNum() != 6 or atom_j.GetAtomicNum() != 6:
                continue
            # Check if either end is close (within 4 bonds) to any indole atom.
            attached = False
            for idx in indole_atoms:
                if dm[i][idx] <= 4 or dm[j][idx] <= 4:
                    attached = True
                    break
            if attached:
                candidate_bonds.append((i,j))
                
    if not candidate_bonds:
        return False, "No attached non‐aromatic vinyl fragment (terpenoid moiety proxy) detected"
    
    # For each candidate we now “grow” the fragment outside the indole core.
    for i, j in candidate_bonds:
        external_atoms = set()
        # If an atom of the double bond is not in the indole, start BFS from it.
        if i not in indole_atoms:
            external_atoms |= bfs_fragment(i, indole_atoms)
        if j not in indole_atoms:
            external_atoms |= bfs_fragment(j, indole_atoms)
        # Count carbon atoms in the external fragment.
        ext_carbons = sum(1 for idx in external_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if ext_carbons >= 5:
            terpene_found = True
            break

    if not terpene_found:
        return False, "No attached non‐aromatic vinyl fragment with sufficient terpenoid fragment size detected"
    
    # ----- 3. Basic size and composition checks -----
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, "Not enough carbon atoms to harbor both tryptophan and a monoterpenoid unit"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a typical monoterpenoid indole alkaloid"
    
    return True, "Contains an indole core and an attached terpene-related fragment of adequate size"

# (Optional) To run tests, you may call the function with SMILES strings, e.g.:
# print(is_monoterpenoid_indole_alkaloid("C/C=C\\1/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]"))