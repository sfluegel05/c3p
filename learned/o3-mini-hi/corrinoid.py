"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: a corrinoid
Definition:
  A corrinoid is a derivative of the corrin nucleus, which contains four reduced or 
  partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one 
  direct carbon-carbon bond linking alpha positions.
  
Heuristics used in this implementation:
  1. Look for a macrocyclic ring (using all rings from the molecule) with:
       • Total ring size between min_size and max_size (we use 13 and 21)
       • Exactly 4 nitrogen atoms.
       • At least 3 and no more than 10 conjugated (double or aromatic) bonds.
  2. Otherwise, search for pyrrole‐like rings using the SMARTS "[nH]1cccc1". 
     Union the atoms from matches that are valid pyrrole rings (having exactly one N)
     and check if the connected union has a total size between min_size and max_size, 
     exactly 4 nitrogens, and between 3 and 10 conjugated bonds.
  3. As a bonus: if cobalt is present and some pyrrole rings were found, we note it 
     (but cobalt alone does not qualify a molecule as a corrinoid).
  
If none of these criteria are met, we conclude that a corrin‐like macrocycle was not detected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# SMARTS pattern to capture a pyrrole-like ring 
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
    
    # User-tunable thresholds
    min_size = 13
    max_size = 21
    min_conj = 3
    max_conj = 10

    # Utility: count conjugated bonds between a set of atom indices.
    def count_conjugated_bonds(atom_indices):
        count = 0
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in atom_indices and a2 in atom_indices:
                # Use the bond's conjugation properties: we count aromatic bonds or explicit double bonds.
                if bond.GetIsAromatic() or bond.GetBondType() == Chem.BondType.DOUBLE:
                    count += 1
        return count

    # --------------------------
    # Heuristic 1: Check every ring from full ring info
    # (this may catch rings that are not in the minimal SSSR set)
    # --------------------------
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        ring_size = len(ring)
        if min_size <= ring_size <= max_size:
            # Count number of nitrogens in the ring.
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 4:
                conj_bonds = count_conjugated_bonds(set(ring))
                if conj_bonds < min_conj:
                    continue  # Not enough conjugation.
                if conj_bonds > max_conj:
                    continue  # Too many conjugated bonds (likely a porphyrin-like scaffold).
                return True, (f"Found a macrocycle (size {ring_size}) with exactly 4 nitrogens "
                              f"and {conj_bonds} conjugated bonds.")
    
    # --------------------------
    # Heuristic 2: Union of pyrrole-like rings 
    # --------------------------
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_smarts)
    valid_pyrrole_sets = []
    for match in pyrrole_matches:
        atoms_in_match = [mol.GetAtomWithIdx(idx) for idx in match]
        if sum(1 for atom in atoms_in_match if atom.GetAtomicNum() == 7) == 1:
            valid_pyrrole_sets.append(set(match))
            
    # We expect 4 pyrrole rings in a corrinoid.
    if valid_pyrrole_sets:
        # Combine all atoms from each pyrrole ring.
        union_atoms = set()
        for ps in valid_pyrrole_sets:
            union_atoms |= ps
        
        union_size = len(union_atoms)
        n_union_nitrogen = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        conj_bonds_union = count_conjugated_bonds(union_atoms)
        
        # Check connectivity of union_atoms.
        union_atoms_list = list(union_atoms)
        visited = set()
        if union_atoms_list:
            stack = [union_atoms_list[0]]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in union_atoms and nbr_idx not in visited:
                        stack.append(nbr_idx)
        fully_connected = (visited == union_atoms)
        
        if fully_connected and (min_size <= union_size <= max_size) and (n_union_nitrogen == 4):
            if conj_bonds_union < min_conj:
                return False, (f"Connected pyrrole macrofragment found (size {union_size}) with 4 nitrogens, "
                               f"but only {conj_bonds_union} conjugated bonds (need at least {min_conj}).")
            if conj_bonds_union > max_conj:
                return False, (f"Connected pyrrole macrofragment found (size {union_size}) with 4 nitrogens, "
                               f"but too many conjugated bonds ({conj_bonds_union} > {max_conj}) – likely not a corrinoid.")
            return True, (f"Found a connected macrofragment from pyrrole rings (size {union_size}) with 4 nitrogens "
                          f"and {conj_bonds_union} conjugated bonds.")
    
    # --------------------------
    # Bonus: if cobalt is present and some pyrrole rings were detected, note extra support.
    # --------------------------
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms and valid_pyrrole_sets:
        return False, ("Cobalt is present and pyrrole-like rings were detected, but macrocycle evidence "
                       "(i.e. overall ring size, connectivity, and conjugation count of 4 nitrogens) remains insufficient.")
    
    return False, "Could not detect a corrin-like macrocycle (neither full ring nor combined pyrrole evidence are conclusive)."


# Example usage:
if __name__ == "__main__":
    # Example: cob(II)yrinic acid c monoamide
    smiles_example = "[H][C@@]12[C@H](CC(O)=O)[C@@](C)(CCC(O)=O)C3=C(C)C4=[N+]5C(=CC6=[N+]7C(=C(C)C8=[N+]([C@]1(C)[C@@](C)(CC(O)=O)[C@@H]8CCC(O)=O)[Co--]57N23)[C@@](C)(CC(N)=O)[C@@H]6CCC(O)=O)C(C)(C)[C@@H]4CCC(O)=O"
    flag, reason = is_corrinoid(smiles_example)
    print("Corrinoid?", flag)
    print("Reason:", reason)