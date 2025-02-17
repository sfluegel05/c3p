"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: A D-hexose is a hexose that has D-configuration at position 5.
In aldoses the stereochemistry (by convention “D”) is determined by the configuration 
of the C5 center (the ring carbon bearing the exocyclic CH2OH substituent) which, 
when correctly assigned by CIP rules, should be “R”.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is defined as a molecule with exactly 6 carbon atoms that contains a sugar ring
    (typically a five-membered furanose or six-membered pyranose ring with one oxygen)
    and the ring carbon bearing the exocyclic CH2OH substituent (C5) should have R configuration.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for the classification decision
    """
    # Parse SMILES string to obtain an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that stereochemistry is assigned (needed for CIP assignment). 
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check that the molecule has exactly 6 carbons.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) != 6:
        return False, "Molecule does not have exactly 6 carbon atoms (not a hexose)"
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_ring = None
    candidate_C5 = None  # This will be the ring carbon with the exocyclic CH2OH substituent.
    
    # Look for a sugar ring: a ring of 5 or 6 atoms with exactly one oxygen.
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygens in the ring.
        oxygens_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(oxygens_in_ring) != 1:
            continue
        # Now, for each ring carbon, check for an attached exocyclic CH2OH group.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Must be carbon.
            if atom.GetAtomicNum() != 6:
                continue
            # Look at all neighbors not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # We expect an exocyclic substituent that is CH2OH: a carbon with 2 hydrogens and one oxygen.
                if nbr.GetAtomicNum() == 6:
                    # Use GetTotalNumHs() for implicit+explicit hydrogens.
                    if nbr.GetTotalNumHs() != 2:
                        continue
                    # Count oxygen neighbors of this substituent.
                    o_count = 0
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 8 and subnbr.GetIdx() != atom.GetIdx():
                            o_count += 1
                    # A CH2OH group should have exactly one oxygen (the hydroxyl).
                    if o_count == 1:
                        candidate_C5 = atom
                        sugar_ring = ring
                        break
            if candidate_C5 is not None:
                break
        if sugar_ring is not None:
            break
    
    if sugar_ring is None or candidate_C5 is None:
        return False, "No appropriate sugar ring with exocyclic CH2OH substituent (C5 candidate) found"
    
    # Check the chirality (CIP code) of the candidate C5.
    if not candidate_C5.HasProp('_CIPCode'):
        return False, "Stereochemistry was not assigned for candidate C5"
    
    cip = candidate_C5.GetProp('_CIPCode')
    if cip != 'R':
        return False, f"C5 has CIP configuration {cip} (expected R for D-hexose)"
    
    return True, "Molecule is a D-hexose: 6 carbons with a sugar ring and the C5 center in R configuration"

# Example usage:
if __name__ == "__main__":
    # This example uses alpha-D-glucose.
    test_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    result, reason = is_D_hexose(test_smiles)
    print(result, reason)