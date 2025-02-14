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
        
    # Get information about the rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    found = False  # Flag to indicate if tetrahydrofuranone ring is found
    
    # Iterate over all rings in the molecule
    for ring in rings:
        # Consider only five-membered rings
        if len(ring) != 5:
            continue
        
        # Count the number of oxygen atoms in the ring
        o_atoms_in_ring = [atom_idx for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8]
        if len(o_atoms_in_ring) != 1:
            continue  # Skip if the ring doesn't have exactly one oxygen atom (not oxolane)
        
        # Check if any ring atom has a double-bonded oxygen attached (oxo substituent)
        has_oxo = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Check for double bond to oxygen outside the ring
                if (bond.GetBondType() == Chem.BondType.DOUBLE and
                    neighbor.GetAtomicNum() == 8 and
                    neighbor.GetIdx() not in ring):
                    has_oxo = True
                    break
            if has_oxo:
                break
        if has_oxo:
            found = True
            break  # Tetrahydrofuranone ring found
    
    if found:
        return True, "Contains an oxolane ring with an oxo substituent"
    else:
        return False, "Does not contain an oxolane ring with an oxo substituent"