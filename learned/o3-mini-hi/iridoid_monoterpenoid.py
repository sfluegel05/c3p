"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid monoterpenoid
Definition: Monoterpenoids biosynthesized from isoprene units, usually consisting of 
a cyclopentane ring fused to a six-membered oxygen heterocycle. Cleavage of a bond 
in the cyclopentane gives rise to the secoiridoids.
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    This function looks for a fused bicyclic system consisting of a cyclopentane ring
    (5-membered) fused to a six-membered ring that contains at least one oxygen atom.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an iridoid monoterpenoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (returns a tuple of atom index tuples for each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Initialize flag and message variables.
    fused_core_found = False
    reason = "No fused cyclopentane and six-membered oxygen heterocycle core found."

    # Iterate over pairs of rings to detect a fused bicyclic system.
    # We require one ring of size 5 and one ring of size 6 that share at least 2 atoms.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = set(ring_info[i])
            ring2 = set(ring_info[j])
            # Identify if one ring is five-membered and the other is six-membered
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                # Determine the fused portion by finding the intersection
                shared_atoms = ring1.intersection(ring2)
                if len(shared_atoms) >= 2:
                    # Identify the six-membered ring from the pair
                    six_ring = ring1 if len(ring1) == 6 else ring2
                    # Check if the six-membered ring contains at least one oxygen atom.
                    oxygen_in_ring = False
                    for idx in six_ring:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetAtomicNum() == 8:
                            oxygen_in_ring = True
                            break
                    if oxygen_in_ring:
                        fused_core_found = True
                        reason = ("Found fused bicyclic core: a cyclopentane ring and a six-membered "
                                  "oxygen heterocycle sharing {} atoms.".format(len(shared_atoms)))
                        break
        if fused_core_found:
            break

    if fused_core_found:
        return True, reason
    else:
        return False, reason