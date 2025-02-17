"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Defined as: Any member of the class of chlorobenzenes carrying three chloro substituents 
at unspecified positions on an aromatic six‐membered ring.
Additional criteria:
  - Reject molecules with multiple fragments (e.g. salts).
  - Reject molecules containing metal atoms.
  - On the candidate benzene ring (six-membered, all carbon, aromatic),
    count chlorine substituents (atomic number 17). All other substituents must be “simple”
    meaning that if the substituent is a single atom it must be one of C, N, O, or S;
    if it is a multi‐atom group then (as in biphenyls) it must be part of another six‐membered
    aromatic (benzene) ring.
"""

from rdkit import Chem

# Define a set of metals to be rejected.
METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
}

# Allowed atomic numbers for a terminal substituent (one‐atom group)
ALLOWED_TERMINAL = {6, 7, 8, 16}

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene.
    A trichlorobenzene is defined as a molecule containing (at least one)
    aromatic six-membered ring (with all atoms being carbon) that
    carries exactly three chlorine substituents; additional substituents
    attached to that ring must be "simple": either a terminal atom (with atomic number
    in ALLOWED_TERMINAL) or be part of another benzene ring [i.e. a six-membered aromatic ring
    formed exclusively of carbon].
    
    Args:
        smiles (str): The SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a trichlorobenzene, False otherwise.
        str: A reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Reject molecules with multiple fragments (salts, counter ions)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (possible salt or counter ion present)."

    # Reject molecules containing metal atoms.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in METALS:
            return False, f"Molecule contains a metal: {atom.GetSymbol()}"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Iterate over each ring in the molecule.
    for ring in atom_rings:
        # Interested only in six-membered rings.
        if len(ring) != 6:
            continue

        # Check that every atom in the ring is aromatic and is a carbon (benzene ring)
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() and atom.GetSymbol() == "C" for atom in ring_atoms):
            continue

        chloro_count = 0
        ring_valid = True  # whether this ring qualifies

        # Loop through each atom in the candidate benzene ring
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check each neighbor of the ring atom that is not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # ignore atoms inside the ring
                # If neighbor is chlorine, count it.
                if nbr.GetAtomicNum() == 17:
                    chloro_count += 1
                else:
                    # Otherwise check if the neighbor is a "simple" substituent.
                    # Case 1: Terminal substituent (degree==1)
                    if nbr.GetDegree() == 1:
                        if nbr.GetAtomicNum() not in ALLOWED_TERMINAL:
                            ring_valid = False
                            break
                    else:
                        # Case 2: Multi-atom substituent.
                        # Allow if the neighbor is part of a six-membered aromatic ring
                        # composed solely of carbon (i.e. another benzene ring; as in biphenyl connections).
                        nbr_rings = nbr.GetOwningMol().GetRingInfo().AtomRings()
                        found_valid_ring = False
                        for subring in nbr_rings:
                            if len(subring) == 6:
                                # Check if this ring (which nbr is part of) is benzene.
                                ring_ok = True
                                for sub_idx in subring:
                                    sub_atom = mol.GetAtomWithIdx(sub_idx)
                                    if not (sub_atom.GetIsAromatic() and sub_atom.GetSymbol() == "C"):
                                        ring_ok = False
                                        break
                                if ring_ok:
                                    found_valid_ring = True
                                    break
                        if not found_valid_ring:
                            ring_valid = False
                            break
            if not ring_valid:
                break

        # Finally, if this ring has exactly three chlorine substituents and passed all allowed tests, we classify it.
        if ring_valid and chloro_count == 3:
            return True, "Found an aromatic six‐membered benzene ring with exactly three chlorine substituents (and all other substituents are simple – terminal atoms or linked benzene rings)."

    return False, "No qualifying aromatic six‐membered benzene ring with exactly three chlorine substituents found."


# Example usage (for testing purposes only – these lines can be removed when integrating into a larger system)
if __name__ == "__main__":
    # A small selection of SMILES strings for testing.
    test_smiles = [
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene (should be True)
        "Clc1ccc(Cl)c(-c2ccccc2)c1Cl",  # 2,3,6-trichlorobiphenyl (should be True)
        "Clc1ccc(Cl)c(Cl)c1",  # 1,2,4-trichlorobenzene (should be True)
        "ClC1=C(Cl)C(OC)=C(Cl)C(=C1OC2=C(Cl)C(OC)=C(Cl)C(=C2OC3=C(Cl)C(OC)=C(Cl)C(=C3O)Cl)Cl)Cl",  # Chlorophellin C (known false positive in previous version)
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")