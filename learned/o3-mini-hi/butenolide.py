"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide – a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.
We look for a five‐membered ring with exactly one ring oxygen and a candidate carbon that bears an exocyclic carbonyl bond.
We further require that the exocyclic oxygen is not involved in any ring (i.e. is truly exocyclic) and that the candidate
either participates in an unsaturated (double) bond with a ring neighbor OR the candidate carbon is sp2-hybridized.
These additional constraints should reduce the false positive rate seen with our earlier attempt.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (gamma-lactone with 2-furanone character or a substituted derivative)
    based on its SMILES string.
    
    The routine checks for a five-membered ring meeting these criteria:
      - Contains exactly one oxygen (the ring-embedded heteroatom).
      - Has at least one carbon (candidate) that bears exactly one exocyclic double bond to an oxygen 
        (the carbonyl) where that oxygen is not a ring member.
      - In addition, the candidate carbon must either be engaged in a double bond with one of its ring-neighbors
        (implying unsaturation) or be sp2-hybridized.
        
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
    ring_atom_sets = ring_info.AtomRings()  # each tuple is a set of indices for one ring
    
    # Loop over all rings looking for a 5-membered ring with exactly one oxygen atom.
    for ring in ring_atom_sets:
        if len(ring) != 5:
            continue
        
        # Count ring oxygen atoms.
        ring_oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if ring_oxygen_count != 1:
            continue
        
        # Now look at each atom in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We expect the candidate to be a carbon.
            if atom.GetAtomicNum() != 6:
                continue
            
            # Look for an exocyclic double bond from this carbon.
            exocyclic_carbonyl_bonds = []
            for bond in atom.GetBonds():
                nb = bond.GetOtherAtom(atom)
                # Only consider bonds going outside the ring:
                if nb.GetIdx() in ring:
                    continue
                # Look for a double bond to oxygen.
                if nb.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    exocyclic_carbonyl_bonds.append(bond)
            
            # We require exactly one exocyclic C=O on this candidate atom.
            if len(exocyclic_carbonyl_bonds) != 1:
                continue
            
            # Check that the exocyclic oxygen is not part of any ring
            oxy = exocyclic_carbonyl_bonds[0].GetOtherAtom(atom)
            in_a_ring = False
            for ring_atoms in ring_atom_sets:
                if oxy.GetIdx() in ring_atoms:
                    in_a_ring = True
                    break
            if in_a_ring:
                continue  # the oxygen must be exocyclic
            
            # Now test for the "furanone character": either
            # (a) the candidate forms a double bond (unsaturation) with at least one ring neighbor, OR
            # (b) the candidate atom is sp2-hybridized.
            unsaturation_in_ring = False
            for nb in atom.GetNeighbors():
                if nb.GetIdx() not in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
                    unsaturation_in_ring = True
                    break
            
            if unsaturation_in_ring or (atom.GetHybridization() == rdchem.HybridizationType.SP2):
                return True, ("Contains a 5-membered gamma-lactone ring with one ring oxygen, an exocyclic carbonyl "
                              "attached to a ring carbon that shows unsaturation or is sp2, consistent with a butenolide")
            # If candidate fails, try the next candidate atom.
        # End loop over atoms in ring.
    # End loop over rings.
    
    return False, "No suitable 5-membered gamma-lactone (butenolide) ring detected"
    
# For testing purposes, you can uncomment the following lines:
# test_smiles = [
#     "[H][C@@]1(CC[C@]2(O)[C@]3([H])CC[C@]4(O)C[C@H](C[C@@H](O)[C@]4(CO)[C@@]3([H])[C@H](O)C[C@]12C)O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C1COC(=O)C1",  # dihydroouabain (True)
#     "O=C1O/C(=C/C2=CC=CC=C2)/C(=C1CC3=CC=CC=C3)CC4=CC=CC=C4",   # Maculalactone C (True)
#     "O=C1C2=C(O)C=3C(=O)C=C(OC)C(C3C=C2C(=O)C4=C1C5=C(C(=O)OC5)C(=C4)OC)=O",  # Paramagnetoquinone B (should be False)
# ]
# for smi in test_smiles:
#     res, reason = is_butenolide(smi)
#     print(smi, res, reason)