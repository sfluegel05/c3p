"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoid – members of the class of tetrahydrochromenochromene compounds.
A rotenoid is defined as a cis‐fused tetrahydrochromeno[3,4-b]chromene skeleton 
(a fused bicyclic system where one ring shows chromenone (benzopyranone) features 
and is fused with another oxygen‐containing ring) and its substituted derivatives.
This implementation uses a heuristic approach based on splitting the fused system 
into (a) a six-membered “chromenone‐like” ring (with exactly one carbonyl attached 
to a ring carbon and a ring oxygen) and (b) a six-membered “ether” ring (with at least 
one oxygen and no carbonyl), plus molecular weight constraints.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule’s molecular weight should be in a range typical for rotenoids (here:
         250–600 Da).
      2. In the ring system the algorithm looks for two fused six-membered rings:
          • One candidate “chromenone” ring that contains at least one oxygen in the ring and 
            exactly one carbon (in the ring) that forms a double bond to an oxygen (the carbonyl).
          • One candidate “ether” ring containing a ring oxygen but no carbonyl groups.
      3. The two candidate rings must be fused – i.e. share at least 2 atoms.
      4. When the two rings are merged the total number of distinct carbonyl groups (from 
         any ring carbon) should equal one.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is determined to be a rotenoid, False otherwise.
        str: A text explanation of the decision.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight (typical rotenoids are between about 250 and 600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (250 <= mol_wt <= 600):
        return False, f"Molecular weight {mol_wt:.1f} out of typical rotenoid range (250–600 Da)"
    
    # Get all ring atom indices (each ring is a tuple of atom indices)
    ring_infos = mol.GetRingInfo().AtomRings()
    
    # Lists to store candidate rings as sets of atom indices
    chromenone_rings = []  # six-membered rings with exactly one C=O (from a ring-carbon) & an oxygen in ring
    ether_rings = []      # six-membered rings with at least one oxygen and NO C=O in ring
    
    # Loop over each ring in the molecule
    for ring in ring_infos:
        if len(ring) != 6:
            continue  # Only process six-membered rings
        
        # Check if there is at least one oxygen in the ring (atom atomic number 8)
        has_ring_oxygen = any(mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring)
        if not has_ring_oxygen:
            continue
        
        # Count how many ring carbons are attached (by a double bond) to an oxygen (carbonyl group)
        carbonyl_atoms = set()
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() != 6:
                continue  # only consider carbon atoms
            # Loop over bonds
            for bond in atom.GetBonds():
                # Check if the bond is a double bond
                if bond.GetBondTypeAsDouble() == 2:
                    # Check if one partner is an oxygen (atomic number 8)
                    other = bond.GetOtherAtom(atom)
                    if other.GetAtomicNum() == 8:
                        carbonyl_atoms.add(i)
                        break  # count only once per carbon
        
        n_carbonyl = len(carbonyl_atoms)
        # Classify the ring:
        if n_carbonyl == 1:
            chromenone_rings.append(set(ring))
        elif n_carbonyl == 0:
            ether_rings.append(set(ring))
        # If there are >1 carbonyls in a ring, it is less likely to be the rotenoid core,
        # so we do not consider it as a candidate.
    
    # Now, look for a pair of candidate rings (one chromenone and one ether) that are fused.
    fused_found = False
    for ring1 in chromenone_rings:
        for ring2 in ether_rings:
            # Fused means at least 2 atoms are shared
            if len(ring1.intersection(ring2)) >= 2:
                # For extra confidence, merge the two rings and count the total carbonyls in the union.
                fused_atoms = ring1.union(ring2)
                overall_carbonyl_atoms = set()
                for i in fused_atoms:
                    atom = mol.GetAtomWithIdx(i)
                    if atom.GetAtomicNum() != 6:
                        continue
                    for bond in atom.GetBonds():
                        if bond.GetBondTypeAsDouble() == 2:
                            other = bond.GetOtherAtom(atom)
                            if other.GetAtomicNum() == 8:
                                overall_carbonyl_atoms.add(i)
                                break
                if len(overall_carbonyl_atoms) == 1:
                    fused_found = True
                    break
        if fused_found:
            break

    if fused_found:
        return True, "Found fused pair: one six‐membered chromenone-like ring (with one carbonyl) fused to a six‐membered ether ring; overall molecular properties are consistent with a rotenoid skeleton."
    else:
        return False, "No characteristic fused chromenone/ether ring pair with a single carbonyl was found; molecule is unlikely to be a rotenoid."

# The entire module is contained in one code block.
# Example test (you may run these to see the reasoning):
if __name__ == "__main__":
    # Example: 13alpha-Hydroxydolineone (a known rotenoid)
    test_smiles = "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C=C6OC=CC6=C5"
    result, reason = is_rotenoid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)