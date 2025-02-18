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
    if it is a multi‐atom group then it must be part of another six‐membered aromatic (benzene) ring.
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
ALLOWED_TERMINAL = {6, 7, 8, 16}  # C, N, O, S

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene.
    
    A trichlorobenzene is defined as a molecule containing at least one
    six-membered aromatic ring (with all atoms being carbon) that carries exactly
    three chlorine substituents. Additionally, any substituent attached to that ring
    must be 'simple': if it is a single atom, it must be one of {C, N, O, S}; if it is
    a multi-atom group it must be part of another benzene ring.
    
    Args:
        smiles (str): A SMILES string of the molecule.
    
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
        # Only interested in six-membered rings.
        if len(ring) != 6:
            continue

        # Check that every atom in the ring is aromatic and is a carbon (benzene ring).
        candidate_ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() and atom.GetSymbol() == "C" for atom in candidate_ring_atoms):
            continue

        chloro_count = 0
        ring_valid = True  # whether this candidate ring qualifies

        # Process each atom in the candidate ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Examine each substituent atom (neighbors not in candidate ring)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # ignore atoms that are part of the candidate ring
                
                # If neighbor is chlorine, count it.
                if nbr.GetAtomicNum() == 17:
                    chloro_count += 1
                else:
                    # For substituents that are not chlorine, they must be "simple"
                    # Case 1: Terminal substituent (one-atom group).
                    if nbr.GetDegree() == 1:
                        if nbr.GetAtomicNum() not in ALLOWED_TERMINAL:
                            ring_valid = False
                            break
                    else:
                        # Case 2: Multi-atom substituent.
                        # Check if the neighbor attaches to a six-membered aromatic ring of carbons
                        # (other than the candidate ring). We use the molecule's ring info and check
                        # only rings that contain the neighbor atom.
                        valid_substituent = False
                        for subring in atom_rings:
                            if nbr.GetIdx() not in subring:
                                continue
                            if set(subring) == set(ring):
                                # This ring is the candidate ring. Skip it.
                                continue
                            if len(subring) == 6:
                                # Verify that all atoms in this subring are aromatic carbons.
                                if all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in subring):
                                    valid_substituent = True
                                    break  # valid benzene ring found for the substituent.
                        if not valid_substituent:
                            ring_valid = False
                            break  # this substituent is not simple
            if not ring_valid:
                break

        # If the candidate ring passed the substituent test and has exactly 3 chlorine atoms, classify as trichlorobenzene.
        if ring_valid and chloro_count == 3:
            return True, ("Found an aromatic six‐membered benzene ring with exactly three chlorine substituents "
                          "and all other substituents are simple (terminal atoms or part of another benzene ring).")
    
    return False, "No qualifying aromatic six‐membered benzene ring with exactly three chlorine substituents found."

# Example usage (for testing purposes only – these lines can be removed when integrating into a larger system)
if __name__ == "__main__":
    test_smiles = [
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene (should be True)
        "Clc1ccc(Cl)c(-c2ccccc2)c1Cl",  # 2,3,6-trichlorobiphenyl (should be True)
        "Clc1ccc(Cl)c(Cl)c1",  # 1,2,4-trichlorobenzene (should be True)
        "ClC1=C(Cl)C(OC)=C(Cl)C(=C1OC2=C(Cl)C(OC)=C(Cl)C(=C2OC3=C(Cl)C(OC)=C(Cl)C(=C3O)Cl)Cl)Cl",  # Expected false positive in the previous version.
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")