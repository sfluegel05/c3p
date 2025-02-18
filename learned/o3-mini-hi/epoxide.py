"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide
Definition: Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
In a valid epoxide, the ring contains exactly one oxygen and two carbons.
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if the molecule represented by the SMILES string is an epoxide.
    An epoxide is defined as a cyclic ether with a 3-membered ring (one oxygen, two carbons).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an epoxide, False otherwise.
        str: Reason for the classification.
    """
    # Convert SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    epoxide_count = 0  # Counter for epoxide rings found.
    # Loop through each ring in the molecule.
    for ring in atom_rings:
        if len(ring) != 3:
            continue  # Only interested in 3-membered rings.
        oxygen_count = 0
        carbon_count = 0
        # Loop over each atom index in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
        # If the 3-membered ring has exactly one oxygen and two carbons, count it as an epoxide.
        if oxygen_count == 1 and carbon_count == 2:
            epoxide_count += 1
            
    if epoxide_count > 0:
        return True, f"Found epoxide ring(s): {epoxide_count} occurrence(s) of a 3-membered cyclic ether"
    else:
        return False, "No epoxide ring (3-membered cyclic ether with one oxygen and two carbons) found"