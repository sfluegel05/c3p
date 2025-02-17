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
  - For SSSR rings, instead of accepting only rings of size 15 or 19, we relax
    the union-of-pyrrole fragments approach to allow a connected macrofragment of
    size 13 to 21 atoms (to catch slight variations).
  - We add a check that within the putative macrocycle there are at least three conjugated (double or aromatic) bonds.
  - Cobalt is used only as extra evidence when there are at least three pyrrole-like fragments.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# SMARTS to capture a pyrrole-like ring (five-membered ring containing a protonated N)
pyrrole_smarts = Chem.MolFromSmarts("[nH]1cccc1")

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    
    Heuristics:
      1. Use default SSSR to look for a macrocycle of size 15 or 19 with exactly 4 nitrogen atoms.
         (If that is found, then we assume a corrinoid core is present.)
      2. Otherwise, search for pyrrole-like rings with the above SMARTS.
         Merge the atoms from all such rings and, if the resulting union is connected,
         has a size between 13 and 21 atoms, and exactly 4 nitrogen atoms, then also check
         that at least three bonds inside this union are either double bonds or aromatic.
      3. As a bonus hint, if cobalt (atomic number 27) is present along with at least 3 pyrrole-like rings,
         that supports (but does not by itself determine) a corrinoid assignment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple with classification (True=corrinoid, False=not) and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --------------------------
    # Heuristic 1: SSSR macrocycle check
    # --------------------------
    sssr = Chem.GetSymmSSSR(mol)
    for ring in sssr:
        ring_atoms = list(ring)
        ring_size = len(ring_atoms)
        if ring_size in (15, 19):
            n_nitrogen = sum(1 for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 4:
                return True, f"Found an SSSR macrocycle of size {ring_size} with exactly 4 nitrogens."

    # --------------------------
    # Heuristic 2: Pyrrole-based union approach
    # --------------------------
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_smarts)
    pyrrole_rings = []
    for match in pyrrole_matches:
        # Double-check that the match (typically 5 atoms) contains exactly one nitrogen.
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        n_N = sum(1 for atom in atoms if atom.GetAtomicNum() == 7)
        if n_N == 1:
            pyrrole_rings.append(set(match))
    
    if len(pyrrole_rings) >= 4:
        union_atoms = set()
        for ring in pyrrole_rings:
            union_atoms |= ring
        # Check connectivity on the union (graph search using bonds only among union atoms).
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
                atom = mol.GetAtomWithIdx(current)
                for nbr in atom.GetNeighbors():
                    idx = nbr.GetIdx()
                    if idx in union_atoms and idx not in visited:
                        stack.append(idx)
            connected = (visited == union_atoms)
        
        union_size = len(union_atoms)
        n_union_nitrogen = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)

        # Check for conjugation within the union: count bonds (between union atoms) that are double or aromatic.
        conjugated_bonds = 0
        # For each bond, check if both endpoints are in union_atoms.
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in union_atoms and a2 in union_atoms:
                bt = bond.GetBondType()
                if bt == Chem.BondType.DOUBLE or bond.GetIsAromatic():
                    conjugated_bonds += 1

        if connected and union_size >= 13 and union_size <= 21 and n_union_nitrogen == 4:
            if conjugated_bonds >= 3:
                return True, (f"Found a connected macrofragment from pyrrole-like rings of size {union_size} "
                              f"with 4 nitrogens and {conjugated_bonds} conjugated bonds.")
            else:
                return False, (f"Found a connected macrofragment from pyrrole-like rings of size {union_size} with 4 nitrogens, "
                               f"but only {conjugated_bonds} conjugated bonds (needed ≥3).")

    # --------------------------
    # Bonus: Cobalt evidence (supportive but not decisive)
    # --------------------------
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms and len(pyrrole_rings) >= 3:
        return False, ("Cobalt is present and some pyrrole-like rings were detected, "
                       "but the macrocycle evidence remains insufficient for a corrinoid assignment.")

    return False, "Could not detect a corrin-like macrocycle (union of pyrrole-like rings with proper size/conjugation) or sufficient evidence."

# Example usage:
if __name__ == "__main__":
    # Example SMILES: cob(II)yrinic acid c monoamide
    smiles_example = "[H][C@@]12[C@H](CC(O)=O)[C@@](C)(CCC(O)=O)C3=C(C)C4=[N+]5C(=CC6=[N+]7C(=C(C)C8=[N+]([C@]1(C)[C@@](C)(CC(O)=O)[C@@H]8CCC(O)=O)[Co--]57N23)[C@@](C)(CC(N)=O)[C@@H]6CCC(O)=O)C(C)(C)[C@@H]4CCC(O)=O"
    flag, reason = is_corrinoid(smiles_example)
    print(flag, reason)