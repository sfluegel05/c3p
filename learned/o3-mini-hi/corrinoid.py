"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: a corrinoid
Definition: A derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings
joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.

This revised implementation uses a two‐step heuristic:
  1. Scan all rings (as computed by RDKit) and look for one macrocycle of size either 15 or 19 that contains exactly 4 nitrogens.
  2. If none is found, then look for 5-membered “pyrrole‐like” rings (with exactly one nitrogen each). Then combine
     these rings into their union and, if this produced a connected cyclic fragment of size 15 or 19 with exactly 4 nitrogens,
     classify as a corrinoid.
In addition, if the molecule contains cobalt and there are at least 3 pyrrole‐like rings, this is used as supportive evidence.
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.

    Heuristics used:
      1. Look for a macrocycle (as given by the SSSR of RDKit) that is either 15 or 19 atoms in size and contains exactly 4 nitrogens.
      2. Otherwise, look for 5-membered rings (“pyrrole‐like rings”) that have exactly one nitrogen.
         Then combine all such rings and, if their union forms a single connected cyclic fragment 
         with a total of 15 or 19 atoms and exactly 4 nitrogen atoms, classify as corrinoid.
      3. As bonus evidence, if cobalt (atomic number 27) is present along with at least 3 pyrrole‐like rings, return True.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a corrinoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # -------------------------------
    # Heuristic 1: Look for a single macrocycle of size 15 or 19 with exactly 4 nitrogens.
    for ring in rings:
        ring_size = len(ring)
        if ring_size in (15, 19):
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 4:
                return True, f"Found a macrocycle of size {ring_size} with exactly 4 ring nitrogens."
    
    # -------------------------------
    # Heuristic 2: Look for pyrrole-like rings: 5-membered rings with exactly one nitrogen.
    pyrrole_rings = []
    for ring in rings:
        if len(ring) == 5:
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 1:
                pyrrole_rings.append(set(ring))
    
    # If we found at least 4 such rings, try to see if they are connected.
    if len(pyrrole_rings) >= 4:
        # Form the union of all atom indices in these rings.
        union_atoms = set()
        for ring in pyrrole_rings:
            union_atoms |= ring

        # Check connectivity for the union: build a simple connectivity check using bonds in mol.
        # (Consider only bonds that connect two atoms in the union.)
        # Use DFS from an arbitrary atom in the union.
        union_atoms_list = list(union_atoms)
        if not union_atoms_list:
            connected = False
        else:
            visited = set()
            stack = [union_atoms_list[0]]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                # Iterate over neighbors of current that are in union_atoms.
                atom = mol.GetAtomWithIdx(current)
                for nbr in atom.GetNeighbors():
                    idx_nbr = nbr.GetIdx()
                    if idx_nbr in union_atoms and idx_nbr not in visited:
                        stack.append(idx_nbr)
            connected = (visited == union_atoms)

        if connected:
            union_size = len(union_atoms)
            n_union_nitrogen = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if union_size in (15, 19) and n_union_nitrogen == 4:
                return True, f"Found a connected macrocycle from pyrrole-like rings of size {union_size} with 4 nitrogens."
    
    # -------------------------------
    # Bonus: Many corrinoids (but not all) include a cobalt ion.
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms and len(pyrrole_rings) >= 3:
        return True, "Cobalt present and several pyrrole-like rings detected, supporting a corrinoid assignment."
    
    return False, "Could not detect a corrin-like macrocycle (15 or 19 atoms with 4 nitrogens) or sufficient connected pyrrole-like rings."

# Example usage (uncomment the following lines to test):
# smiles_example = "[H][C@]12[C@H](CC(O)=O)[C@@](C)(CCC(O)=O)C3=C(C)C4=[N+]5C(=CC6=[N+]7C(=C(C)C8=[N+]([C@]1(C)[C@@](C)(CC(O)=O)[C@@H]8CCC(O)=O)[Co--]57N23)[C@@](C)(CC(N)=O)[C@@H]6CCC(O)=O)C(C)(C)[C@@H]4CCC(O)=O"  # cob(II)yrinic acid c monoamide
# flag, reason = is_corrinoid(smiles_example)
# print(flag, reason)