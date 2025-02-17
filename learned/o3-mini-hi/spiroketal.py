"""
Classifies: CHEBI:72600 spiroketal
"""
#!/usr/bin/env python
"""
Classifies: spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
"""

from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal features a ketal carbon that has two oxygen substituents and is the only atom shared between two rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a spiroketal, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of rings, each is a tuple of atom indices
    
    # If no rings, it cannot be spiroketal
    if not rings:
        return False, "No rings detected in the molecule"
    
    # Loop over all atoms looking for a spiro carbon candidate
    for atom in mol.GetAtoms():
        # Check only carbons as candidate ketal centers.
        if atom.GetAtomicNum() != 6:
            continue
        
        atom_idx = atom.GetIdx()
        
        # Collect all rings that contain this atom
        rings_for_atom = [set(ring) for ring in rings if atom_idx in ring]
        
        # Need at least two rings to be a spiro center
        if len(rings_for_atom) < 2:
            continue
        
        # Check every pair of rings; if they share only this atom, we may have a spiro center.
        spiro_found = False
        for i in range(len(rings_for_atom)):
            for j in range(i + 1, len(rings_for_atom)):
                # The spiro condition: the intersection of the two rings is exactly this atom.
                if rings_for_atom[i].intersection(rings_for_atom[j]) == {atom_idx}:
                    spiro_found = True
                    break
            if spiro_found:
                break
                
        # If a potential spiro center was found, check if it is a ketal carbon.
        if spiro_found:
            # Count the oxygen substituents bonded to this candidate carbon.
            oxygen_bonds = 0
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    oxygen_bonds += 1
            if oxygen_bonds >= 2:
                return True, f"Found spiroketal center at atom index {atom_idx} with {oxygen_bonds} oxygen substituents"
            else:
                # Even though the atom is a spiro center, it does not show the ketal oxygen pattern.
                return False, f"Found spiro carbon at atom index {atom_idx}, but it has only {oxygen_bonds} oxygen substituents (need >=2)"
    
    return False, "No spiroketal pattern found in the molecule"

# Example usage (these lines can be commented out or removed in production):
if __name__ == "__main__":
    test_smiles = "[C@]12([C@@H]([C@](O[C@]3(O1)[C@H]([C@H](O[C@]3(CC2=O)C)C)C(=O)O)([C@H](C)/C=C(\\C)/C=C/C(=C/4\\C(NCC4=O)=O)/O)[H])C)[H]"  # nocamycin E example
    result, reason = is_spiroketal(test_smiles)
    print(result, reason)