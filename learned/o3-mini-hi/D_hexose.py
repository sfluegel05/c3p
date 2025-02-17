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
    A D-hexose is defined here as a molecule with 6 carbon atoms that contains a sugar ring
    (typically a five-membered furanose or six-membered pyranose ring with one oxygen)
    and that has the chiral center corresponding to C5 (the ring carbon bearing an exocyclic CH2OH group)
    with R configuration.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for the classification decision
    """
    # Parse the SMILES string to obtain an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry (required for CIP assignment)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for hexose: must contain exactly 6 carbons.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) != 6:
        return False, "Molecule does not have exactly 6 carbon atoms (not a hexose)"
    
    # Next, try to detect a sugar ring: a ring that is either 5 or 6 members in size
    # (furanose or pyranose) and that contains exactly one oxygen.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_ring = None
    candidate_C5 = None  # This will be our candidate for the 'C5' stereocenter.
    
    # Loop over rings and look for one with 5 or 6 atoms and exactly one oxygen.
    for ring in rings:
        if len(ring) not in (5, 6):
            continue  # not a typical sugar ring size
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        oxygens_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(oxygens_in_ring) != 1:
            continue  # sugar rings in hexoses have one ring oxygen
        # Now search for a ring carbon that has an exocyclic CH2OH group.
        # In most cyclic hexoses the CH2OH group is attached to C5.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # we need a carbon
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() != 6:
                    continue
                # Check if the neighbor is CH2OH. We expect it to be methylene (CH2) with an -OH.
                h_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                # Count oxygen neighbors (apart from the connection back to the ring)
                o_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetIdx() != atom.GetIdx())
                if h_count == 2 and o_count >= 1:
                    candidate_C5 = atom
                    sugar_ring = ring
                    break
            if candidate_C5 is not None:
                break
        if sugar_ring is not None:
            break
    
    if sugar_ring is None or candidate_C5 is None:
        return False, "No appropriate sugar ring with exocyclic CH2OH substituent (C5 candidate) found"
    
    # Check the chirality (CIP) of the candidate C5.
    # For D-hexoses the C5 stereocenter (the one carrying the -CH2OH substituent) should have the R configuration.
    if not candidate_C5.HasProp('_CIPCode'):
        return False, "Stereochemistry was not assigned for candidate C5"
    cip = candidate_C5.GetProp('_CIPCode')
    if cip != 'R':
        return False, f"C5 has CIP configuration {cip} (expected R for D-hexose)"
    
    return True, "Molecule is a D-hexose: 6 carbons with a sugar ring and C5 in R configuration"

# Example usage:
if __name__ == "__main__":
    test_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"  # example: alpha-D-glucose
    result, reason = is_D_hexose(test_smiles)
    print(result, reason)