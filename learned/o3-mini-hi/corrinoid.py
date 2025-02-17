"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: a corrinoid
Definition:
  A corrinoid is a derivative of the corrin nucleus, which contains four reduced or 
  partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one 
  direct carbon-carbon bond linking alpha positions.
  
Improvements compared to a previous attempt:
  - Uses RDKit's GetSymmSSSR to help capture ring information.
  - Searches for a macrocycle (either from SSSR or via a union of pyrrole-like rings)
    that is either 15 or 19 atoms in size with exactly 4 nitrogens.
  - Uses a SMARTS for pyrrole-like rings instead of only ring size.
  - Applies cobalt presence as bonus evidence only if clear pyrrole fragments are found.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# SMARTS for a pyrrole-like ring: five-membered ring with one nitrogen (with at least one bonded hydrogen)
pyrrole_smarts = Chem.MolFromSmarts("[nH]1cccc1")

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    
    Heuristics:
      1. Check for a macrocycle from ring data using GetSymmSSSR that has size 15 or 19 and exactly 4 nitrogens.
      2. If not found, search for pyrrole-like rings using a SMARTS filter. 
         Merge the atom indices from such rings and check whether the union is a connected cycle
         of size 15 or 19 having exactly 4 nitrogen atoms.
      3. Only as bonus, if cobalt (atomic number 27) is present AND at least 3 pyrrole-like rings were detected,
         then that evidence supports a corrinoid assignment. (But it is not used alone.)
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a corrinoid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Precompute ring information using SSSR (the default symmetry SSSR)
    sssr = Chem.GetSymmSSSR(mol)
    rings = [tuple(r) for r in sssr]  # each ring is a tuple of atom indices

    # Heuristic 1:
    for ring in rings:
        ring_size = len(ring)
        if ring_size in (15, 19):
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 4:
                return True, f"Found an SSSR macrocycle of size {ring_size} with exactly 4 nitrogens."

    # Heuristic 2: Use pyrrole SMARTS.
    # Find all substructure matches of a pyrrole-like ring.
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_smarts)
    pyrrole_rings = []
    for match in pyrrole_matches:
        # Although the SMARTS itself finds 5 atoms, double-check that there is exactly one nitrogen.
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        if sum(1 for atom in atoms if atom.GetAtomicNum() == 7) == 1:
            pyrrole_rings.append(set(match))
    
    # If there are at least 4 pyrrole-like rings, try to see if they connect into a convex fragment.
    if len(pyrrole_rings) >= 4:
        union_atoms = set()
        for ring in pyrrole_rings:
            union_atoms |= ring
        
        # Check connectivity of the union (using only bonds between atoms in union_atoms).
        union_atoms_list = list(union_atoms)
        if union_atoms_list:
            visited = set()
            stack = [union_atoms_list[0]]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                atom = mol.GetAtomWithIdx(current)
                for nbr in atom.GetNeighbors():
                    idx_nbr = nbr.GetIdx()
                    if idx_nbr in union_atoms and idx_nbr not in visited:
                        stack.append(idx_nbr)
            connected = (visited == union_atoms)
        else:
            connected = False
        
        union_size = len(union_atoms)
        n_union_nitrogen = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if connected and union_size in (15, 19) and n_union_nitrogen == 4:
            return True, f"Found a connected macrocycle from pyrrole-like rings of size {union_size} with 4 nitrogens."
    
    # Bonus: Only consider cobalt as supportive evidence if there is at least modest support from pyrrole rings.
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms and len(pyrrole_rings) >= 3:
        # But if neither Heuristic 1 nor 2 picked it up, then cobalt alone is not enough.
        return False, "Cobalt present and some pyrrole-like rings detected, but insufficient macrocycle evidence."
    
    return False, "Could not detect a corrin-like macrocycle (15 or 19 atoms with 4 nitrogens) or sufficient connected pyrrole-like rings."

# Example usage:
if __name__ == "__main__":
    # Example: cob(II)yrinic acid c monoamide
    smiles_example = "[H][C@@]12[C@H](CC(O)=O)[C@@](C)(CCC(O)=O)C3=C(C)C4=[N+]5C(=CC6=[N+]7C(=C(C)C8=[N+]([C@]1(C)[C@@](C)(CC(O)=O)[C@@H]8CCC(O)=O)[Co--]57N23)[C@@](C)(CC(N)=O)[C@@H]6CCC(O)=O)C(C)(C)[C@@H]4CCC(O)=O"
    flag, reason = is_corrinoid(smiles_example)
    print(flag, reason)