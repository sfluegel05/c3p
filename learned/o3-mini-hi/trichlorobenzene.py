"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Defined as: any member of the class of chlorobenzenes carrying three chloro substituents 
at unspecified positions on an aromatic six‐membered ring (and that, aside from the chlorine, 
all other substituents on such ring are attached via carbon or oxygen).
Additional criteria: 
  - Reject molecules with multiple fragments (e.g. salts).
  - Reject molecules containing metal atoms.
"""

from rdkit import Chem

# Define a set of metal element symbols that we want to reject.
METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
}

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene.
    A trichlorobenzene is defined as a molecule containing a six-membered aromatic ring 
    that has exactly three chlorine atoms attached as substituents and that, aside from 
    chlorine, any other substituents attached directly to the ring come only from carbon or oxygen.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a trichlorobenzene, False otherwise.
        str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules with multiple fragments (e.g. salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (possible salt or counter ion present)"

    # Reject molecules containing metal atoms.
    for atom in mol.GetAtoms():
        # Instead of atom.GetIsMetal(), we use the atom's symbol and our METALS set.
        if atom.GetSymbol() in METALS:
            return False, f"Molecule contains a metal: {atom.GetSymbol()}"

    # Obtain ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Loop over each ring in the molecule
    for ring in atom_rings:
        # We are interested only in 6-membered rings.
        if len(ring) != 6:
            continue

        # Check if all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        chloro_count = 0
        ring_ok = True  # flag to ensure no disallowed substituents on the ring

        # For each atom in the ring, examine neighboring atoms which are not part of the ring.
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # ignore atoms that are part of the ring
                atomic_num = nbr.GetAtomicNum()
                # Count chlorine substituents (atomic number 17).
                if atomic_num == 17:
                    chloro_count += 1
                else:
                    # For non-chlorine substituents, allow only if they are attached via carbon (6) or oxygen (8).
                    if atomic_num not in (6, 8):
                        ring_ok = False
                        break
            if not ring_ok:
                break

        if not ring_ok:
            # Skip rings that include disallowed substituents.
            continue

        # If a ring has exactly three chlorines then it qualifies as a trichlorobenzene.
        if chloro_count == 3:
            return True, "Found an aromatic six-membered ring with exactly three chloro substituents (all other substituents via C or O)."

    return False, "No qualifying aromatic six-membered ring with exactly three chloro substituents found."


# Example usage (these prints can be removed or commented out when integrating into a larger system):
if __name__ == "__main__":
    # Test SMILES strings: you can add more examples as needed.
    test_smiles = [
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene; valid
        "Clc1ccc(Cl)c(-c2ccccc2)c1Cl",  # 2,3,6-trichlorobiphenyl; not a simple benzene ring (has an extra phenyl group)
        "C1(=CC=C(C(=C1Cl)CC(O)=O)Cl)Cl"  # chlorfenac; has multiple fragments or disallowed substituents
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")