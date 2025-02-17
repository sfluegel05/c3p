"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: a corrinoid
Definition:
  A corrinoid is a derivative of the corrin nucleus, which contains four reduced or 
  partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one 
  direct carbon-carbon bond linking alpha positions.
  
Heuristics:
  1. Look for an SSSR ring with size between 13 and 21 atoms that has exactly 4 nitrogen atoms.
  2. Otherwise, search for pyrrole-like rings (SMARTS: [nH]1cccc1), union their atoms, and if the resulting
     connected subgraph has a total number of atoms between 13 and 21 with exactly 4 nitrogens and at least 3
     conjugated (double or aromatic) bonds, the molecule is classified as a corrinoid.
  3. Additionally, if cobalt (atomic number 27) is present together with at least 3 pyrrole-like rings, this lends
     extra support (though not conclusive on its own) to the corrinoid assignment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# SMARTS pattern to capture a pyrrole-like ring (typically a 5-membered ring with a protonated nitrogen)
pyrrole_smarts = Chem.MolFromSmarts("[nH]1cccc1")

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the boolean indicates whether the molecule is classified as a corrinoid,
                     and the string explains the reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --------------------------
    # Heuristic 1: SSSR macrocycle check (relaxed to any ring of size 13 to 21 with 4 nitrogens)
    # --------------------------
    sssr_rings = Chem.GetSymmSSSR(mol)
    for ring in sssr_rings:
        ring_atoms = list(ring)
        ring_size = len(ring_atoms)
        # Accept any ring whose size is between 13 and 21
        if 13 <= ring_size <= 21:
            n_nitrogen = sum(1 for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 4:
                # Also count conjugated (double or aromatic) bonds in the ring
                conjugated_bonds = 0
                for bond in mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx()
                    a2 = bond.GetEndAtomIdx()
                    if a1 in ring_atoms and a2 in ring_atoms:
                        if bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetIsAromatic():
                            conjugated_bonds += 1
                if conjugated_bonds >= 3:
                    return True, (f"Found an SSSR macrocycle of size {ring_size} with exactly 4 nitrogens "
                                  f"and {conjugated_bonds} conjugated bonds.")
    
    # --------------------------
    # Heuristic 2: Union of pyrrole-like rings
    # --------------------------
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_smarts)
    pyrrole_rings = []
    for match in pyrrole_matches:
        # Ensure that each 5-atom match indeed has exactly one nitrogen.
        atoms_in_match = [mol.GetAtomWithIdx(idx) for idx in match]
        if sum(1 for atom in atoms_in_match if atom.GetAtomicNum() == 7) == 1:
            pyrrole_rings.append(set(match))
    
    # Only proceed if we found at least 4 candidate pyrrole rings (we expect 4 pyrrole rings in a corrinoid core)
    if len(pyrrole_rings) >= 4:
        union_atoms = set()
        for ring in pyrrole_rings:
            union_atoms |= ring
        
        # Check connectivity of atoms in union_atoms
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
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in union_atoms and nbr_idx not in visited:
                        stack.append(nbr_idx)
            connected = (visited == union_atoms)
        else:
            connected = False
        
        union_size = len(union_atoms)
        n_union_nitrogen = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        
        # Count conjugated bonds within the union atoms (double bonds or aromatic bonds)
        conjugated_bonds = 0
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in union_atoms and a2 in union_atoms:
                if bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetIsAromatic():
                    conjugated_bonds += 1
        
        if connected and (13 <= union_size <= 21) and n_union_nitrogen == 4:
            if conjugated_bonds >= 3:
                return True, (f"Found a connected macrofragment from pyrrole-like rings of size {union_size} "
                              f"with 4 nitrogens and {conjugated_bonds} conjugated bonds.")
            else:
                return False, (f"Found a connected macrofragment from pyrrole-like rings of size {union_size} with 4 nitrogens, "
                               f"but only {conjugated_bonds} conjugated bonds (need at least 3).")
    
    # --------------------------
    # Bonus: Check if cobalt is present and at least 3 pyrrole-like rings were detected.
    # (This fact supports, but does not fulfill, corrinoid classification.)
    # --------------------------
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms and len(pyrrole_rings) >= 3:
        return False, ("Cobalt is present and some pyrrole-like rings were detected, "
                       "but macrocycle evidence (size, connectivity, conjugation, 4 nitrogens) remains insufficient.")
    
    return False, "Could not detect a corrin-like macrocycle (neither SSSR nor combined pyrrole evidence are conclusive)."


# Example usage:
if __name__ == "__main__":
    # Example SMILES: cob(II)yrinic acid c monoamide
    smiles_example = "[H][C@@]12[C@H](CC(O)=O)[C@@](C)(CCC(O)=O)C3=C(C)C4=[N+]5C(=CC6=[N+]7C(=C(C)C8=[N+]([C@]1(C)[C@@](C)(CC(O)=O)[C@@H]8CCC(O)=O)[Co--]57N23)[C@@](C)(CC(N)=O)[C@@H]6CCC(O)=O)C(C)(C)[C@@H]4CCC(O)=O"
    flag, reason = is_corrinoid(smiles_example)
    print("Corrinoid?", flag)
    print("Reason:", reason)