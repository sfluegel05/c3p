"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoids – members of the class of tetrahydrochromenochromene compounds.
A rotenoid is defined as having a cis‐fused tetrahydrochromeno[3,4‐b]chromene skeleton.
This skeletal motif comprises two six‐membered rings:
  • a “chromenone‐like” ring that carries exactly one exocyclic carbonyl (attached to a ring carbon) 
    and at least one ring oxygen (but not necessarily exactly one, to allow minor substituent effects),
    and is not fully aromatic;
  • an “ether” ring that has at least one oxygen and no exocyclic carbonyl, and also is not fully aromatic.
The two rings must be fused by sharing exactly 2 adjacent atoms (which must be carbons).
Additional substituents (e.g. glycosides) may be present without affecting the rotenoid core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    
    Heuristic criteria (improved):
      1. Molecule must have a minimum molecular weight (>=200 Da).
      2. Look for all six‐membered rings that contain at least one oxygen and are not fully aromatic.
         For every such ring, count how many exocyclic carbonyls are present on carbons in that ring.
         (An exocyclic carbonyl is defined as a double–bond from a ring carbon to an oxygen not in the ring.)
         Then define:
           • Candidate chromenone‐like ring: at least one ring oxygen with exactly one exocyclic carbonyl.
           • Candidate ether ring: at least one ring oxygen with zero exocyclic carbonyl.
      3. Look for a fused pair (one candidate from each list) that share exactly two atoms.
         Moreover, require that in each ring the shared atoms appear consecutively in the (cyclic) ring order,
         and that both shared atoms are carbons.
      4. For the union of the two rings, ensure that there is exactly one exocyclic carbonyl overall.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple with boolean classification and textual explanation.
                   If the rotenoid criteria are met, returns True with an explanation;
                   otherwise returns False with a reason.
    """
    
    # Helper: Check if in a cyclic ordered list the given set of two indices appear consecutively
    def shared_atoms_adjacent(ordered_ring, shared_set):
        n = len(ordered_ring)
        for i in range(n):
            pair = {ordered_ring[i], ordered_ring[(i+1) % n]}
            if pair == shared_set:
                return True
        return False

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Minimal molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} is too low for rotenoids"
    
    # Retrieve ring info (RDKit returns tuples of atom indices for every ring)
    ring_infos = mol.GetRingInfo().AtomRings()
    
    chromenone_candidates = []  # Each entry: (ordered_ring, set(ring))
    ether_candidates = []       # Each entry: (ordered_ring, set(ring))
    
    # Loop over all rings and consider only six-membered rings
    for ring in ring_infos:
        if len(ring) != 6:
            continue
        # Get atom objects for the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Must contain at least one oxygen
        if not any(atom.GetAtomicNum() == 8 for atom in ring_atoms):
            continue
        # Skip if the ring is fully aromatic (we expect partial saturation)
        if all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Count how many oxygens are in the ring (only within the ring indices)
        ring_oxygen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        
        # Count the number of exocyclic carbonyls attached to ring carbons.
        exo_carbonyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbons for exocyclic carbonyl search
            if atom.GetAtomicNum() != 6:
                continue
            # Check bonds of this carbon for a double bond to an oxygen not in the ring.
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8 and (neighbor.GetIdx() not in ring):
                        exo_carbonyl_count += 1
                        break  # one per ring atom max
        
        # Classify the ring as candidate chromenone-like if it has at least one oxygen and exactly one exocyclic carbonyl;
        # candidate ether ring if it has at least one oxygen and no exocyclic carbonyl.
        if exo_carbonyl_count == 1:
            chromenone_candidates.append( (ring, set(ring)) )
        elif exo_carbonyl_count == 0:
            ether_candidates.append( (ring, set(ring)) )
        # Rings with >1 exocyclic carbonyl are not considered.
    
    if not chromenone_candidates:
        return False, "No candidate chromenone-like ring (6-membered with ≥1 ring oxygen and one exocyclic carbonyl) found"
    if not ether_candidates:
        return False, "No candidate ether ring (6-membered with ≥1 ring oxygen and no exocyclic carbonyl) found"
    
    # Look for a fused pair of candidate rings (one chromenone-like and one ether) that share exactly 2 atoms.
    for chrom_ring, chrom_set in chromenone_candidates:
        for ether_ring, ether_set in ether_candidates:
            shared_atoms = chrom_set.intersection(ether_set)
            if len(shared_atoms) != 2:
                continue
            # Check that both shared atoms are carbons.
            if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in shared_atoms):
                continue
            # Check that in both rings the two shared atoms appear consecutively (cyclic order).
            if not (shared_atoms_adjacent(chrom_ring, shared_atoms) and shared_atoms_adjacent(ether_ring, shared_atoms)):
                continue
            
            # Now examine the union (fused core) to ensure overall only one exocyclic carbonyl is present.
            fused_union = chrom_set.union(ether_set)
            union_exo_carbonyls = 0
            for idx in fused_union:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for bond in atom.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8 and (neighbor.GetIdx() not in fused_union):
                            union_exo_carbonyls += 1
                            break
            if union_exo_carbons == 1 or union_exo_carbonyls == 1:
                reason_msg = ("Found fused pair: one six‐membered chromenone‐like ring (with one exocyclic carbonyl and at least one ring oxygen) "
                              "cis‐fused (sharing exactly 2 adjacent carbons) to a six‐membered ether ring; molecular properties are consistent "
                              "with a rotenoid skeleton.")
                return True, reason_msg
            
    return False, "No characteristic cis‐fused chromenone/ether ring pair (sharing exactly 2 adjacent carbons with the expected exocyclic carbonyl pattern) was found; molecule is unlikely to be a rotenoid."

# Example test code (can be removed or executed)
if __name__ == "__main__":
    # Test one of the true positive examples: 13alpha-Hydroxydolineone
    test_smiles = "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C=C6OC=CC6=C5"
    result, reason = is_rotenoid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)