"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is a lactam in which the amide bond is contained within a four-membered ring,
    which includes the amide nitrogen and the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate if beta-lactam ring is found
    beta_lactam_found = False

    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if ring is a 4-membered ring
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_types = [atom.GetAtomicNum() for atom in atoms_in_ring]

            # Check if the ring contains exactly one nitrogen atom and three carbon atoms
            num_nitrogens = atom_types.count(7)  # Atomic number 7 for nitrogen
            num_carbons = atom_types.count(6)    # Atomic number 6 for carbon

            if num_nitrogens == 1 and num_carbons == 3:
                # Identify the nitrogen atom in the ring
                nitrogen_atom = next(atom for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
                # Check if nitrogen has a bond to a carbonyl carbon within the ring
                for bond in nitrogen_atom.GetBonds():
                    neighbor = bond.GetOtherAtom(nitrogen_atom)
                    if neighbor.GetIdx() in ring and neighbor.GetAtomicNum() == 6:
                        # Check if this carbon is a carbonyl carbon (double bonded to oxygen)
                        is_carbonyl = False
                        for bond2 in neighbor.GetBonds():
                            if bond2.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                oxy_neighbor = bond2.GetOtherAtom(neighbor)
                                if oxy_neighbor.GetAtomicNum() == 8:  # Oxygen
                                    if oxy_neighbor.GetIdx() in ring:
                                        is_carbonyl = True
                                        break
                        if is_carbonyl:
                            beta_lactam_found = True
                            break
                if beta_lactam_found:
                    break

    if beta_lactam_found:
        return True, "Molecule contains a beta-lactam ring"
    else:
        return False, "No beta-lactam ring found"