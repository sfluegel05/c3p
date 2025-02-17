"""
Classifies: CHEBI:72600 spiroketal
"""
#!/usr/bin/env python
"""
Classifies: spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
The ketal carbon in a spiroketal is expected to be sp³-hybridized, have exactly four substituents,
with exactly two oxygen atoms and two carbon atoms attached, and it must lie at the junction
(where a pair of rings have no other atoms in common).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    
    A spiroketal features a ketal (sp³) carbon that is the only common atom bridging two rings.
    Here we require that candidate ketal carbons
     - Are sp³-hybridized carbons.
     - Belong to at least two rings.
     - Have at least one pair of rings whose only intersection is that carbon.
     - Have exactly four bonds, of which exactly two are to oxygen atoms and two to carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the molecule"
    
    # Loop over all atoms searching for a candidate spiroketal carbon.
    for atom in mol.GetAtoms():
        # Consider only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check if the carbon is sp3 hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        
        atom_idx = atom.GetIdx()
        # Get the rings that include this atom.
        rings_for_atom = [set(ring) for ring in rings if atom_idx in ring]
        if len(rings_for_atom) < 2:
            continue  # Needs to be part of at least two rings.
        
        # Check for at least one pair of rings that share only this atom.
        is_spiro_center = False
        for i in range(len(rings_for_atom)):
            for j in range(i + 1, len(rings_for_atom)):
                if rings_for_atom[i].intersection(rings_for_atom[j]) == {atom_idx}:
                    is_spiro_center = True
                    break
            if is_spiro_center:
                break
        
        if not is_spiro_center:
            continue
        
        # Now verify that this carbon shows the ketal pattern:
        # It should be tetravalent with exactly 4 bonds,
        # exactly 2 bonds to oxygen atoms, and the remainder to carbon atoms.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 4:
            continue  # Not tetravalent; likely not a ketal center.
        
        oxygen_count = 0
        carbon_count = 0
        for nbr in neighbors:
            atomic_num = nbr.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
            else:
                # If other atom types are involved, then it might not be a standard ketal.
                pass
        
        if oxygen_count == 2 and carbon_count == 2:
            return True, f"Found spiroketal center at atom index {atom_idx} with {oxygen_count} oxygen and {carbon_count} carbon substituents"
        # Otherwise continue and check other candidate atoms.
    
    return False, "No spiroketal pattern found in the molecule"

# Example usage:
if __name__ == "__main__":
    # Test with the nocamycin E SMILES given in the problem statement.
    test_smiles = "[C@]12([C@@H]([C@](O[C@]3(O1)[C@H]([C@H](O[C@]3(CC2=O)C)C)C(=O)O)([C@H](C)/C=C(\\C)/C=C/C(=C/4\\C(NCC4=O)=O)/O)[H])C)[H]"
    result, reason = is_spiroketal(test_smiles)
    print(result, reason)