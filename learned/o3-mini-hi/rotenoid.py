"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoids – members of the class of tetrahydrochromenochromene compounds.
A rotenoid is defined as having a cis‐fused tetrahydrochromeno[3,4‐b]chromene skeleton.
This consists of two six‐membered rings fused by sharing exactly two atoms.
The first ring (“chromenone‐like”) must contain at least one ring oxygen and exactly one carbonyl
(double‐bonded oxygen attached to a ring carbon, where that oxygen is exocyclic).
The second ring (“ether” ring) must contain at least one oxygen and no carbonyl groups.
Additional substituents (e.g. glycosides) may be present.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    
    Heuristic criteria (improved):
      1. The molecule must have a molecular weight >= 200 Da.
      2. Identify all six-membered rings that contain at least one oxygen.
         For each such ring:
           a. If exactly one of the ring carbons has a double bond to an oxygen (exocyclic), 
              mark it as a candidate chromenone-like ring.
           b. If no ring carbon has a double bond to oxygen, mark it as a candidate ether ring.
      3. Search for a pair (one from each candidate list) that are fused by exactly 2 atoms.
         Moreover, require the shared atoms are both carbons.
      4. Optionally, check that the union of the two rings does not show additional carbonyls beyond
         that single carbonyl in the chromenone ring.
         
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple with boolean classification and a reason.
                   If criteria are met, returns True with an explanation;
                   otherwise, returns False with the reason.
    """
    
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for minimal molecular weight (>=200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} is too low for rotenoids"
    
    # Retrieve ring information from the molecule
    ring_infos = mol.GetRingInfo().AtomRings()
    
    # Lists to store candidate rings (as sets of atom indices)
    chromenone_candidates = []  # should have exactly one carbonyl (double bond from a ring C) 
    ether_candidates = []       # no carbonyl on any ring carbon
    
    # Loop through each ring; only consider six-membered rings with at least one oxygen.
    for ring in ring_infos:
        if len(ring) != 6:
            continue
            
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Proceed only if the ring contains at least one oxygen atom (atomic num 8)
        if not any(atom.GetAtomicNum() == 8 for atom in ring_atoms):
            continue
        
        carbonyl_count = 0  # count ring carbons with a double bond to oxygen (exocyclic oxygen)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # For every bond of this ring carbon check if there is a double bond to oxygen.
            for bond in atom.GetBonds():
                # Check for a double bond.
                if bond.GetBondTypeAsDouble() == 2:
                    # Identify the neighboring atom.
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        # Check if the oxygen is exocyclic; i.e. not part of the ring.
                        if neighbor.GetIdx() not in ring:
                            carbonyl_count += 1
                            break  # count each ring carbon only once
        
        # Classify ring based on carbonyl_count.
        if carbonyl_count == 1:
            chromenone_candidates.append(set(ring))
        elif carbonyl_count == 0:
            ether_candidates.append(set(ring))
        # Rings with >1 carbonyl (or other counts) are not considered as part of the core.
    
    if not chromenone_candidates:
        return False, "No candidate chromenone-like ring (6-membered with one exocyclic carbonyl and ring oxygen) found"
    
    if not ether_candidates:
        return False, "No candidate ether ring (6-membered with oxygen and no in-ring carbonyl) found"
    
    # Search for a fused pair of rings that share exactly 2 atoms.
    fused_pair_found = False
    reason_msg = ""
    for chrom_ring in chromenone_candidates:
        for ether_ring in ether_candidates:
            shared_atoms = chrom_ring.intersection(ether_ring)
            if len(shared_atoms) == 2:
                # Additional check: ensure both shared atoms are carbons.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in shared_atoms):
                    # Merge the rings and check total carbonyls in the fused union.
                    fused_union = chrom_ring.union(ether_ring)
                    union_carbonyls = 0
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
                    # We expect one carbonyl (from the chromenone) in the fused core.
                    if union_carbonyls == 1:
                        fused_pair_found = True
                        reason_msg = ("Found fused pair: one six‐membered chromenone-like ring (with one exocyclic carbonyl "
                                      "and ring oxygen) cis-fused (sharing exactly 2 carbon atoms) to a six‐membered ether ring; "
                                      "molecular properties are consistent with a rotenoid skeleton.")
                        break
        if fused_pair_found:
            break
    
    if fused_pair_found:
        return True, reason_msg
    else:
        return False, "No characteristic cis-fused chromenone/ether ring pair (sharing exactly 2 carbons) was found; molecule is unlikely to be a rotenoid."

# Example test code
if __name__ == "__main__":
    # Known rotenoid: 13alpha-Hydroxydolineone
    test_smiles = "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C=C6OC=CC6=C5"
    result, reason = is_rotenoid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)