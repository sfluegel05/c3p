"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide – a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.

Our approach:
  - Parse the molecule and retrieve all rings.
  - Loop over every 5‐membered ring that contains exactly one ring oxygen.
  - For that ring, check every carbon neighbor of the ring oxygen.
  - For each such candidate carbon determine if it bears an exocyclic carbonyl group (i.e. a double bond to an oxygen that is not in the ring).
  - If such a motif is found we classify the structure as butenolide.
  
Note: We have relaxed the original requirements (e.g. internal unsaturation of the candidate carbon)
because many butenolide derivatives (e.g. dihydroouabain, neriifolin) do not strictly satisfy that constraint.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (a gamma-lactone with a 2-furanone skeleton or its substituted/dihydro derivatives)
    based on its SMILES string.

    The function proceeds as follows:
      - Parses the molecule and obtains all rings.
      - For each 5-membered ring with exactly one oxygen, it examines each ring carbon directly bonded to the ring oxygen.
      - It then looks for an exocyclic double-bond to oxygen (i.e. a keto carbonyl) on that candidate carbon.
      - If such a pattern is found the molecule is classified as a butenolide.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule meets the butenolide criteria; False otherwise.
        str: A message with the reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get information about all rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over every ring looking for a 5-membered ring with exactly one oxygen.
    for ring in rings:
        if len(ring) != 5:
            continue  # skip rings that are not 5-membered

        # Count oxygen atoms in the ring.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # require exactly one oxygen in the ring

        # Identify the lactone (ring) oxygen.
        lactone_oxygen_idx = ring_oxygens[0]
        lactone_oxygen = mol.GetAtomWithIdx(lactone_oxygen_idx)
        
        # Look at neighbors of the lactone oxygen that lie in the ring (candidate carbons).
        for neighbor in lactone_oxygen.GetNeighbors():
            if neighbor.GetIdx() not in ring:
                continue
            if neighbor.GetAtomicNum() != 6:
                continue  # we want a carbon neighbor
            
            candidate_carbon = neighbor
            # Look for an exocyclic double bond on this candidate carbon.
            exocyclic_carbonyl_found = False
            for bond in candidate_carbon.GetBonds():
                # Only consider bonds that go outside the ring.
                other_atom = bond.GetOtherAtom(candidate_carbon)
                if other_atom.GetIdx() in ring:
                    continue
                # Check if the bond is a double bond and the other atom is oxygen.
                if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    exocyclic_carbonyl_found = True
                    break
            if exocyclic_carbonyl_found:
                return True, ("Contains a 5‐membered gamma‐lactone (butenolide) ring with one ring oxygen "
                              "and a candidate ring carbon bearing an exocyclic carbonyl group, "
                              "consistent with a 2‐furanone-type motif.")
    
    return False, "No suitable 5‐membered gamma‐lactone (butenolide) ring with the proper 2-furanone motif detected"
    
# Example usage:
# test_smiles = [
#   "[H][C@@]1(CC[C@]2(O)[C@]3([H])CC[C@]4(O)C[C@H](C[C@@H](O)[C@]4(CO)[C@@]3([H])[C@H](O)C[C@]12C)O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C1COC(=O)C1",  # dihydroouabain
#   "O=C1O/C(=C/C2=CC=CC=C2)/C(=C1CC3=CC=CC=C3)CC4=CC=CC=C4",   # Maculalactone C
# ]
# for sm in test_smiles:
#     result, reason = is_butenolide(sm)
#     print(sm, result, reason)