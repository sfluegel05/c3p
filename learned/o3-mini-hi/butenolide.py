"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide – a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.
We look for a five‐membered ring with exactly one ring oxygen in which at least one carbon atom:
  • Bears an exocyclic carbonyl group (a double bond to an O outside the ring), and
  • Either is involved in a double bond with an adjacent ring atom (i.e. unsaturation in the ring)
    or is sp2‐hybridized (allowing for dihydro derivatives).
This is one attempt to improve upon earlier versions that yielded high false positive rates.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    
    A butenolide is defined here as a gamma-lactone with a 2-furanone core or a substituted derivative.
    For our purposes we require that the molecule contains a five-membered ring that:
      - Contains exactly one oxygen (as the ring heteroatom),
      - Has at least one carbon atom of the ring that bears an exocyclic carbonyl group (i.e. a C=O 
        bond with the oxygen located outside the ring),
      - And this candidate carbon is either directly involved in a double bond with a neighboring ring
        atom (i.e. the lactone ring shows unsaturation) or is sp2-hybridized (which may be the case in
        dihydro substituted derivatives).
        
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a butenolide, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    
    # Loop over all rings in the molecule
    for ring in ring_info.AtomRings():
        # Only consider 5-membered rings
        if len(ring) != 5:
            continue
        
        # Count number of oxygen atoms in the ring.
        ring_oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        # Typical 2-furanone rings have exactly one ring oxygen.
        if ring_oxygen_count != 1:
            continue
        
        # For each atom in the ring, check if it is a candidate lactone carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We expect a carbon (atomic number 6) for the carbonyl-bearing atom.
            if atom.GetAtomicNum() != 6:
                continue
            candidate_found = False
            # Check all bonds from this atom that lead OUTSIDE the ring.
            for bond in atom.GetBonds():
                nb = bond.GetOtherAtom(atom)
                if nb.GetIdx() in ring:
                    continue  # only interested in bonds going out of the ring
                # Check for exocyclic double-bond to an oxygen.
                if nb.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    candidate_found = True
                    break
            if not candidate_found:
                continue  # move to next atom if no exocyclic carbonyl is attached
            
            # At this point, we have a candidate: a ring carbon with an exocyclic C=O.
            # Now check for “furanone character”. We require that either:
            # (a) at least one bond within the ring (adjacent to this candidate) is a double bond, or
            # (b) the candidate carbon is sp2-hybridized.
            unsaturation_in_ring = False
            # Get the list of neighbor indices (in the ring only) for the candidate:
            ring_neighbor_idxs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in ring]
            # Check bonds between candidate and each ring neighbor:
            for nb_idx in ring_neighbor_idxs:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb_idx)
                if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    unsaturation_in_ring = True
                    break
            # Also consider the possibility that even if no adjacent bond is double,
            # the candidate atom itself is sp2-hybridized (using RDKit's perception).
            if unsaturation_in_ring:
                return True, ("Contains a 5-membered gamma-lactone ring with an exocyclic carbonyl and "
                              "unsaturation in the ring (2-furanone core)")
            elif atom.GetHybridization() == rdchem.HybridizationType.SP2:
                return True, ("Contains a 5-membered gamma-lactone ring with an exocyclic carbonyl and "
                              "the candidate carbon is sp2 (suggesting a dihydro-2-furanone derivative)")
            # Else: candidate did not pass the extra unsaturation/hybridization check.
        # End for each atom in the current ring.
    # End for each ring.
    
    return False, "No suitable 5-membered gamma-lactone (2-furanone or substituted derivative) ring detected"

# For testing purposes (you can uncomment the following lines):
# test_smiles = [
#     "O=C1O/C(=C/C2=CC=CC=C2)/C(=C1CC3=CC=CC=C3)CC4=CC=CC=C4",   # Maculalactone C (should be True)
#     "[H][C@@]1(CC[C@]2(O)[C@]3([H])CC[C@]4(O)C[C@H](C[C@@H](O)[C@]4(CO)[C@@]3([H])[C@H](O)C[C@]12C)O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C1COC(=O)C1",  # dihydroouabain (should be True)
#     "C1CC2=C(CC1)C(N(C2=O)C3=C(C=C(C(=C3)OC(C)C#C)Cl)F)=O",      # Example false positive structure
# ]
# for smi in test_smiles:
#     res, reason = is_butenolide(smi)
#     print(smi, res, reason)