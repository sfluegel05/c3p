"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: A D-hexose is a hexose that has D-configuration at position 5.
For our purposes a “pure” D-hexose must have exactly 6 carbon atoms, 6 oxygen atoms,
contain a sugar ring (a 5- or 6-membered ring containing exactly one oxygen) and have the ring carbon
bearing a CH2OH substituent (i.e. with a –CH2OH exocyclic group) showing the CIP configuration 'R'.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined here as a molecule that:
      - has exactly 6 carbon atoms and 6 oxygen atoms (as expected for a pure hexose)
      - contains a sugar ring (either a 5-membered or 6-membered ring that contains exactly 1 oxygen atom)
      - in that ring there is one (and only one) ring-carbon that carries an exocyclic CH2OH substituent.
      - the ring-carbon bearing the CH2OH must have CIP stereochemistry R (this is the C5 center).
      
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is classified as a D-hexose, False otherwise.
        str: A message giving the reason for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force assignment of stereochemistry.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Count carbons and oxygens.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(carbons) != 6:
        return False, "Molecule does not have exactly 6 carbon atoms (not a hexose)"
    if len(oxygens) != 6:
        return False, "Molecule does not have exactly 6 oxygen atoms (likely modified hexose)"
    
    # Retrieve ring information and try to identify a sugar ring.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_ring = None
    candidate_C5 = None  # This will be the ring carbon having the CH2OH substituent.
    
    # Look for a ring of size 5 or 6 that has exactly one oxygen atom.
    # (Sugar rings are usually either furanoses (5 members) or pyranoses (6 members)
    for ring in rings:
        if len(ring) not in (5,6):
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_oxygens = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # not typical for a sugar ring
        # Now search the ring atoms (which must be carbons) for an exocyclic substituent that matches CH2OH.
        local_candidate = None
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            if ring_atom.GetAtomicNum() != 6:
                continue  # we expect a carbon
            # Examine neighbors not in the ring.
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # We expect the CH2OH group: a carbon atom (not in the ring) that is a methylene group
                # and that is attached to exactly one oxygen (the hydroxyl) aside from the bond to the ring.
                if nbr.GetAtomicNum() == 6:
                    # Check that the candidate carbon has the right number of hydrogens.
                    # (A –CH2OH group should ideally have 2 hydrogens.)
                    if nbr.GetTotalNumHs() != 2:
                        continue
                    # Count oxygen neighbors (excluding the ring_atom to which it is attached).
                    o_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetIdx() != ring_atom.GetIdx()]
                    if len(o_neighbors) != 1:
                        continue
                    # If we already found a candidate on this ring, then adding another would be ambiguous.
                    if local_candidate is not None:
                        # More than one candidate on the same ring: ambiguous.
                        local_candidate = None
                        break
                    local_candidate = ring_atom
            if local_candidate is not None:
                candidate_C5 = local_candidate
                sugar_ring = ring
                break
        if sugar_ring is not None and candidate_C5 is not None:
            break

    if sugar_ring is None or candidate_C5 is None:
        return False, "No appropriate sugar ring with exocyclic CH2OH substituent (C5 candidate) found"
    
    # Check that the candidate C5 atom has been assigned stereochemistry.
    if not candidate_C5.HasProp('_CIPCode'):
        return False, "Stereochemistry was not assigned for the candidate C5 atom"
    cip = candidate_C5.GetProp('_CIPCode')
    if cip != 'R':
        return False, f"C5 has CIP configuration {cip} (expected R for D-hexose)"
    
    return True, "Molecule is a D-hexose: 6 carbons with a sugar ring and the C5 center in R configuration"

# Example usage:
if __name__ == "__main__":
    # Several SMILES can be tested. Here is a classic case: alpha-D-galactose.
    test_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    result, reason = is_D_hexose(test_smiles)
    print(result, reason)