"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoids – members of the class of tetrahydrochromenochromene compounds.
A rotenoid is defined as having a cis‐fused tetrahydrochromeno[3,4‐b]chromene skeleton.
This skeletal motif comprises two six‐membered rings:
  • a “chromenone‐like” ring having exactly one ring oxygen and exactly one exocyclic carbonyl (double‐bonded oxygen attached to a ring carbon but not in the ring);
  • an “ether” ring that has at least one oxygen and no exocyclic carbonyl.
The two rings must be fused by sharing exactly 2 atoms; moreover these two shared atoms must be adjacent (i.e. the rings are cis‐fused).
Additional substituents (e.g. glycosides) may be present without affecting the rotenoid core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    
    Improved heuristic criteria:
      1. Molecule must have a minimum molecular weight (>=200 Da).
      2. Search for all six‐membered rings that include at least one oxygen.
         For each such ring (in the order given by RDKit):
           a. Count how many ring atoms are oxygens.
           b. For each ring carbon, count if it has an exocyclic carbonyl (a double bond to an oxygen that is not in the ring).
         Then classify:
           • If a ring has exactly 1 oxygen in the ring and exactly 1 exocyclic carbonyl, it is a candidate “chromenone‐like” ring.
           • If a ring has at least 1 oxygen in the ring but no exocyclic carbonyl, it is a candidate “ether” ring.
      3. Look for a pair (one candidate from each list) that share exactly two atoms.
         Moreover, require that these two shared atoms are adjacent 
         (i.e. appear consecutively in the ring order, with wrap‐around taken into account) in both rings 
         and that both shared atoms are carbons.
      4. Additionally, for the union of the two rings, ensure that only one exocyclic carbonyl is present.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple with boolean classification and a textual explanation.
                   If criteria are met, returns True with an explanation;
                   otherwise, returns False with the reason.
    """

    # Helper function: check if two atoms in a ring (given as an ordered tuple) are consecutive (cyclically)
    def shared_atoms_adjacent(ring_order, shared_set):
        n = len(ring_order)
        # check each consecutive pair including wrap-around (i.e. element n-1 and 0)
        for i in range(n):
            pair = {ring_order[i], ring_order[(i+1) % n]}
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

    # Retrieve ring info (RDKit returns ordered tuples of atom indices)
    ring_infos = mol.GetRingInfo().AtomRings()
    
    # Lists for candidate rings
    chromenone_candidates = []  # (ring_order as tuple, set(ring))
    ether_candidates = []       # (ring_order as tuple, set(ring))
    
    # Loop over all rings; consider only six-membered rings with at least one oxygen.
    for ring in ring_infos:
        if len(ring) != 6:
            continue
        # Ensure ring contains at least one oxygen atom
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not any(atom.GetAtomicNum() == 8 for atom in ring_atoms):
            continue

        # Count ring oxygens
        ring_oxygen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        
        # Count exocyclic carbonyls (for ring carbons only; a bond of type DOUBLE between a C in ring and an O not in ring).
        exo_carbonyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # only consider carbons
                continue
            # Examine bonds for a double-bond to an oxygen
            for bond in atom.GetBonds():
                # RDKit represents bond order using GetBondTypeAsDouble (==2.0 for double bond)
                if bond.GetBondTypeAsDouble() == 2:
                    neighbor = bond.GetOtherAtom(atom)
                    # Check that neighbor is oxygen and not part of the ring
                    if neighbor.GetAtomicNum() == 8 and (neighbor.GetIdx() not in ring):
                        exo_carbonyl_count += 1
                        break  # count each ring atom at most once

        # Classify ring based on its features.
        if exo_carbonyl_count == 1 and ring_oxygen_count == 1:
            chromenone_candidates.append( (ring, set(ring)) )
        elif exo_carbonyl_count == 0 and ring_oxygen_count >= 1:
            ether_candidates.append( (ring, set(ring)) )
        # Do not consider rings that do not match either pattern

    if not chromenone_candidates:
        return False, "No candidate chromenone-like ring (6-membered with one ring oxygen and one exocyclic carbonyl) found"
    if not ether_candidates:
        return False, "No candidate ether ring (6-membered with oxygen and no exocyclic carbonyl) found"

    # Search for fused pair: one candidate chromenone ring and one candidate ether ring sharing exactly 2 atoms.
    for chrom_order, chrom_set in chromenone_candidates:
        for ether_order, ether_set in ether_candidates:
            shared_atoms = chrom_set.intersection(ether_set)
            if len(shared_atoms) != 2:
                continue
            # Check that both shared atoms are carbons.
            if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in shared_atoms):
                continue
            # Check if in both rings the shared atoms appear as adjacent (cyclic order).
            if not (shared_atoms_adjacent(chrom_order, shared_atoms) and shared_atoms_adjacent(ether_order, shared_atoms)):
                continue
            # Form the fused union (set of atom indices in either ring)
            fused_union = chrom_set.union(ether_set)
            union_carbonyls = 0
            # For every carbon in the fused union check for an exocyclic carbonyl bond.
            for idx in fused_union:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for bond in atom.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8 and (neighbor.GetIdx() not in fused_union):
                            union_carbonyls += 1
                            break
            # Expect exactly one exocyclic carbonyl (coming from the chromenone ring) in the fused core.
            if union_carbonyls == 1:
                reason_msg = ("Found fused pair: one six‐membered chromenone-like ring (with one exocyclic carbonyl and one ring oxygen) "
                              "cis-fused (sharing exactly 2 adjacent carbons) to a six‐membered ether ring; molecular properties are "
                              "consistent with a rotenoid skeleton.")
                return True, reason_msg

    return False, "No characteristic cis-fused chromenone/ether ring pair (sharing exactly 2 adjacent carbons with the expected exocyclic carbonyl pattern) was found; molecule is unlikely to be a rotenoid."

# Example test code (this section can be removed or executed as desired)
if __name__ == "__main__":
    # Example of a known rotenoid: 13alpha-Hydroxydolineone
    test_smiles = "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C=C6OC=CC6=C5"
    result, reason = is_rotenoid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)