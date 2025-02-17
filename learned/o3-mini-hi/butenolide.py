"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide – a gamma-lactone consisting of a 2-furanone skeleton and its substituted/dihydro derivatives.

This improved version builds on the previous approach but introduces extra filters:
  - After identifying 5-membered rings with exactly one oxygen, it computes the number
    of double bonds within the ring.
  - For each ring carbon neighboring the ring oxygen, it looks for an exocyclic carbonyl group.
  - If the candidate carbon is sp2 and the ring shows no internal double bond (besides the exocyclic C=O),
    the match is less likely to reflect a genuine 2-furanone motif.
  - These additional checks help to filter out false positives.

Note that this approach is heuristic and may not perfectly classify every borderline case.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (a gamma-lactone with a 2-furanone skeleton or its
    substituted/dihydro derivatives) based on its SMILES string.
    
    The function performs the following steps:
      - Parses the molecule.
      - Retrieves all rings and loops over each 5-membered ring containing exactly one oxygen.
      - For each such ring, examines ring carbons that are directly bonded to that oxygen.
      - Checks if a candidate carbon features an exocyclic double bond to oxygen (a carbonyl group).
      - Additionally, it computes the number of double bonds internal to the ring and verifies
        that if the candidate carbon is sp2 hybrids then the ring’s unsaturation is consistent
        with an unsaturated (2-furanone) system; otherwise, it may be a dihydro derivative.
    
    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        bool: True if the molecule meets the butenolide criteria; False otherwise.
        str: Explanation message for the classification decision.
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get information about all rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over every ring; we want 5-membered rings only.
    for ring in rings:
        if len(ring) != 5:
            continue  # not a 5-membered ring
        
        # Count oxygen atoms in the ring.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # require exactly one oxygen atom in the ring
        
        lactone_oxygen_idx = ring_oxygens[0]
        lactone_oxygen = mol.GetAtomWithIdx(lactone_oxygen_idx)
        
        # Count internal (ring-only) double bonds. We only count each bond once.
        ring_double_bonds = 0
        for i in range(len(ring)):
            a1 = ring[i]
            for j in range(i+1, len(ring)):
                a2 = ring[j]
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    ring_double_bonds += 1

        # For each neighbor of the lactone oxygen within the ring, check for an exocyclic C=O
        for neighbor in lactone_oxygen.GetNeighbors():
            if neighbor.GetIdx() not in ring:
                continue  # must be inside the ring
            if neighbor.GetAtomicNum() != 6:
                continue  # candidate should be carbon

            candidate_carbon = neighbor

            # Look for an exocyclic bond (bond to an atom outside of the ring)
            exocyclic_carbonyl_found = False
            for bond in candidate_carbon.GetBonds():
                # Only consider bonds that lead outside the ring
                other_atom = bond.GetOtherAtom(candidate_carbon)
                if other_atom.GetIdx() in ring:
                    continue
                # Check: bond must be double and the other atom should be oxygen
                if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    exocyclic_carbonyl_found = True
                    break
            
            if exocyclic_carbonyl_found:
                # Additional check based on candidate carbon hybridization and ring unsaturation.
                # For a classical 2-furanone, the candidate carbon is sp2 and the ring should have an internal double bond.
                # Allow dihydro derivatives (candidate may be sp3) if no additional double bond is found.
                if candidate_carbon.GetHybridization() == rdchem.HybridizationType.SP2:
                    if ring_double_bonds < 1:
                        # sp2 candidate in a fully saturated ring (apart from the exocyclic bond) is unexpected.
                        continue
                # Accept the candidate match.
                return True, ("Contains a 5‐membered gamma‐lactone (butenolide) ring with one ring oxygen and "
                              "a ring carbon bearing an exocyclic carbonyl group, and with ring unsaturation "
                              f"({ring_double_bonds} internal double bond{'s' if ring_double_bonds != 1 else ''}) "
                              "consistent with a 2‐furanone or dihydro derivative.")
    
    # If no candidate motif is found
    return False, "No suitable 5‐membered gamma‐lactone (butenolide) motif detected"

# Example usage:
# test_smiles = [
#   "[H][C@@]1(CC[C@]2(O)[C@]3([H])CC[C@]4(O)C[C@H](C[C@@H](O)[C@]4(CO)[C@@]3([H])[C@H](O)C[C@]12C)O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C1COC(=O)C1",  # dihydroouabain
#   "O=C1O/C(=C/C2=CC=CC=C2)/C(=C1CC3=CC=CC=C3)CC4=CC=CC=C4",   # Maculalactone C
#   # ... add more test SMILES as needed ...
# ]
# for sm in test_smiles:
#     result, reason = is_butenolide(sm)
#     print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")