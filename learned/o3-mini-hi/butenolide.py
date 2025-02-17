"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide – a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.

Our approach is to (1) find a 5‐membered ring that contains exactly one oxygen atom,
(2) require that the ring oxygen is directly bonded to a carbon that bears an exocyclic C=O (the lactone carbonyl)
    where the oxygen in the carbonyl is not itself in any ring,
(3) further demand that the candidate carbon is sp2‐hybridized and is also engaged in at least one double bond to 
    one of its ring neighbors (i.e. it shows internal unsaturation), in line with the 2–furanone motif.
    
If a ring meeting these criteria is found we classify the structure as a butenolide.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (a gamma-lactone with a 2-furanone skeleton)
    based on its SMILES string.
    
    The function works as follows:
      - Parses the molecule and retrieves all rings.
      - For each five-membered ring with exactly one ring oxygen:
          * Identify the lone ring oxygen.
          * Look at each ring carbon that is directly bonded to this oxygen.
          * For each such candidate carbon, search for an exocyclic double bond to oxygen.
          * Ensure that:
               - There is exactly one exocyclic C=O on that carbon,
               - That the exocyclic oxygen is not a member of any ring,
               - The candidate carbon is sp2-hybridized, and
               - The candidate carbon forms a double bond with at least one ring neighbor.
      - If such a motif is found, the molecule is classified as butenolide.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule meets the butenolide criteria; False otherwise.
        str: A message with the reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    ring_atom_sets = ring_info.AtomRings()  # each tuple is a ring (indices of atoms)
    
    # Loop over all rings looking for a five-membered ring with exactly one oxygen.
    for ring in ring_atom_sets:
        if len(ring) != 5:
            continue
        
        # Count the number of oxygen atoms within the ring.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # skip rings that do not have exactly one ring oxygen
        
        # For the typical 2-furanone motif, the lactone oxygen (in the ring) should be adjacent 
        # to the carbon that bears the exocyclic carbonyl.
        lactone_oxygen_idx = ring_oxygens[0]
        lactone_oxygen = mol.GetAtomWithIdx(lactone_oxygen_idx)
        
        # Get all ring neighbors of the ring oxygen (neighbors that are part of this ring).
        oxygen_ring_neighbors = [nbr for nbr in lactone_oxygen.GetNeighbors() if nbr.GetIdx() in ring]
        
        # Loop over these neighbors that are carbons: candidates for the lactone carbonyl location.
        for candidate in oxygen_ring_neighbors:
            if candidate.GetAtomicNum() != 6:
                continue
            
            # Look for exocyclic double bonds from this candidate carbon.
            exo_carbonyl_bonds = []
            for bond in candidate.GetBonds():
                # Consider only bonds leading outside the ring.
                other = bond.GetOtherAtom(candidate)
                if other.GetIdx() in ring:
                    continue
                # Look for a bond that is a C=O double bond.
                if other.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    exo_carbonyl_bonds.append(bond)
            
            if len(exo_carbonyl_bonds) != 1:
                # We require exactly one exocyclic carbonyl on this candidate.
                continue
            
            # Check that the exocyclic oxygen is not part of any ring.
            exo_oxygen = exo_carbonyl_bonds[0].GetOtherAtom(candidate)
            in_a_ring = False
            for ring_atoms in ring_atom_sets:
                if exo_oxygen.GetIdx() in ring_atoms:
                    in_a_ring = True
                    break
            if in_a_ring:
                continue  # the oxygen must be exocyclic
            
            # Verify that the candidate carbon is sp2 hybridized.
            if candidate.GetHybridization() != rdchem.HybridizationType.SP2:
                continue
            
            # Additionally, check that the candidate carbon forms at least one double bond with a ring neighbor.
            unsaturation_in_ring = False
            for nbr in candidate.GetNeighbors():
                if nbr.GetIdx() not in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(candidate.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    unsaturation_in_ring = True
                    break
            if not unsaturation_in_ring:
                continue
            
            # If we reach here, we found a 5-membered ring with one ring oxygen such that
            # one of its carbons (directly attached to the ring oxygen) bears exactly one exocyclic C=O,
            # is sp2-hybridized and shows at least one internal double bond.
            return True, ("Contains a 5-membered gamma-lactone ring (butenolide) where the lactone oxygen "
                          "is directly bonded to a sp2-hybridized carbon bearing an exocyclic carbonyl and internal unsaturation, "
                          "consistent with a 2-furanone motif.")
    
    return False, "No suitable 5-membered gamma-lactone (butenolide) ring with the proper 2-furanone motif detected"

# For testing purposes, one might run:
# test_smiles = [
#     "[H][C@@]1(CC[C@]2(O)[C@]3([H])CC[C@]4(O)C[C@H](C[C@@H](O)[C@]4(CO)[C@@]3([H])[C@H](O)C[C@]12C)O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C1COC(=O)C1",  # dihydroouabain (True)
#     "O=C1O/C(=C/C2=CC=CC=C2)/C(=C1CC3=CC=CC=C3)CC4=CC=CC=C4",   # Maculalactone C (True)
#     "O=C1C2=C(O)C=3C(=O)C=C(OC)C(C3C=C2C(=O)C4=C1C5=C(C(=O)OC5)C(=C4)OC)=O",  # Paramagnetoquinone B (False)
# ]
# for sm in test_smiles:
#     res, msg = is_butenolide(sm)
#     print(sm, res, msg)