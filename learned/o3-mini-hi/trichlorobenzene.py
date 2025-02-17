"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene 
Defined as: Any member of the class of chlorobenzenes carrying three chloro substituents 
at unspecified positions on a six‐membered aromatic benzene ring. Additional criteria:
  - Reject molecules with multiple fragments (e.g. salts).
  - Reject molecules containing metal atoms.
  - On the candidate benzene ring (six-membered, all carbon, aromatic) count chlorine substituents.
    For each substituent that is not chlorine the “simplicity” is checked as follows:
       • If the substituent (i.e. the fragment attached via a bond to the ring) consists of a single heavy atom,
         it must be one of C, N, O, or S.
       • If the substituent has 2 or 3 heavy atoms, then we allow common small groups: 
             methoxy (C–O), carboxyl ([C(=O)O]) or carboxamide ([C(=O)N]).
       • Alternatively, if the entire substituent is itself an aromatic benzene ring (six-membered, all C),
         it is taken as “simple.”
"""

from rdkit import Chem

# Set of metal symbols to reject.
METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
}

# Allowed terminal atom atomic numbers (if the substituent consists of one atom only).
ALLOWED_TERMINAL = {6, 7, 8, 16}  # corresponding to C, N, O, S

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene as defined above.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule qualifies as a trichlorobenzene, False otherwise.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string."
    
    # Reject molecules with multiple fragments:
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (possible salt or counter ion present)."
    
    # Reject molecules containing metals:
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in METALS:
            return False, f"Molecule contains a metal: {atom.GetSymbol()}"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper: given a starting atom index outside a candidate ring, traverse (without crossing back into the ring)
    # to collect the substituent fragment atoms.
    def get_substituent_fragment(start_idx, ring_set):
        stack = [start_idx]
        visited = set()
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue  # do not cross back into the candidate ring
                if nbr.GetIdx() not in visited:
                    stack.append(nbr.GetIdx())
        return visited  # set of atom indices in the substituent fragment
    
    # Helper: decide if a substituent fragment (given as a set of atom indices) is "simple".
    def is_simple_fragment(fragment_idxs):
        # Get the heavy atom symbols (ignore hydrogens since they are implicit)
        atoms = [mol.GetAtomWithIdx(idx) for idx in fragment_idxs]
        heavy = [atom for atom in atoms if atom.GetAtomicNum() > 1]
        count = len(heavy)
        # If only one heavy atom, allow if it is in ALLOWED_TERMINAL.
        if count == 1:
            if heavy[0].GetAtomicNum() in ALLOWED_TERMINAL:
                return True
            else:
                return False
        # If the fragment is exactly two heavy atoms, allow if it corresponds to a methoxy group: C-O.
        if count == 2:
            nums = sorted([atom.GetAtomicNum() for atom in heavy])
            if nums == [6, 8]:
                return True
        # If the fragment is exactly three heavy atoms, allow if it corresponds to a carboxyl ([C(=O)O])
        # or a carboxamide ([C(=O)N]). We check by sorting atomic numbers.
        if count == 3:
            nums = sorted([atom.GetAtomicNum() for atom in heavy])
            if nums == [6, 8, 8] or nums == [6, 7, 8]:
                return True
        # Lastly, allow if the fragment is itself an aromatic benzene ring (six heavy atoms all being aromatic carbons).
        if count == 6:
            # Check if every atom in the fragment is aromatic and a carbon.
            if all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in fragment_idxs):
                return True
        return False

    # Now iterate over candidate rings.
    for ring in atom_rings:
        # Only consider six-membered rings.
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is aromatic and is a carbon.
        candidate_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() and atom.GetSymbol() == "C" for atom in candidate_atoms):
            continue

        ring_set = set(ring)
        chloro_count = 0   # count chlorine substituents directly attached to the ring
        all_substituents_ok = True
        
        # To avoid double counting substituents (if attached to two ring atoms), keep track of visited bonds.
        checked_bonds = set()
        # For each atom of the ring, check neighbors that are not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                bond_id = tuple(sorted((idx, nbr_idx)))
                if nbr_idx in ring_set or bond_id in checked_bonds:
                    continue
                checked_bonds.add(bond_id)
                if nbr.GetAtomicNum() == 17:  # chlorine
                    chloro_count += 1
                else:
                    frag = get_substituent_fragment(nbr_idx, ring_set)
                    if not is_simple_fragment(frag):
                        all_substituents_ok = False
                        break
            if not all_substituents_ok:
                break
        
        # We require exactly three Cl substituents on the candidate benzene ring.
        if all_substituents_ok and chloro_count == 3:
            return True, ("Found an aromatic six‐membered benzene ring (all C) with exactly three chlorine substituents "
                          "and all other substituents are simple (either a single atom, a small group corresponding to methoxy, "
                          "carboxyl, carboxamide, or another benzene ring).")
    
    return False, "No qualifying aromatic six‐membered benzene ring with exactly three chlorine substituents and allowed substituents found."

# The module can be tested by running a few examples.
if __name__ == "__main__":
    test_smiles = [
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene should be True.
        "Clc1ccc(Cl)c(-c2ccccc2)c1Cl",  # 2,3,6-trichlorobiphenyl should be True.
        "C1=C(C(=C(C(=C1Cl)OC)C(O)=O)Cl)Cl",  # tricamba (should be True by our criteria even though it has extra groups)
        "ClC=1C(C=2C(Cl)=CC(Cl)=C(Cl)C2)=CC(Cl)=C(Cl)C1Cl"  # one of the known false positives – likely False.
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")