"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is defined as 'Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.'
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    found = False  # Flag to indicate if tetrahydrofuranone ring is found

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Iterate over all rings in the molecule
    for ring in rings:
        # Consider only five-membered rings
        if len(ring) != 5:
            continue
        
        # Check if ring contains at least one oxygen atom (oxolane ring)
        oxygen_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(oxygen_in_ring) < 1:
            continue  # Not an oxolane ring
        
        # For each atom in the ring, check if it has a double bond to an oxygen atom
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:  # Neighbor is oxygen
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Confirm that the oxygen is not part of the ring (to avoid counting ring oxygens)
                        if not ring_info.IsAtomInRingOfSize(neighbor.GetIdx(), 5):
                            found = True
                            break
                        else:
                            # Also consider cases where the double-bonded oxygen is in the ring (e.g., lactones)
                            found = True
                            break
            if found:
                break
        if found:
            break

    if found:
        return True, "Contains an oxolane ring with an oxo substituent (tetrahydrofuranone)"
    else:
        return False, "Does not contain an oxolane ring with an oxo substituent"