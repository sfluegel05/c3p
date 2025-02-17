"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoid – members of the class of tetrahydrochromenochromene compounds.
A rotenoid is defined as a cis‐fused tetrahydrochromeno[3,4-b]chromene skeleton (i.e. a fused bicyclic system
where one ring shows chromenone (benzopyranone)‐like features and is fused with another oxygen-containing ring)
and its substituted derivatives.
Note: This implementation uses a heuristic approach based on ring‐info in the molecule.
"""

from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    This heuristic calls a molecule “rotenoid” if it contains a fused bicyclic system made
    of at least two six-membered rings that display key features of a chromenone ring (an ether ring
    with a carbonyl function) fused to another oxygen-containing ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is determined to be a rotenoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    # List to hold candidate six-membered rings with chromenone features.
    candidate_rings = []
    
    # Loop over rings and look for six-membered rings that contain:
    # 1. At least one oxygen atom (present as a heteroatom in the ring).
    # 2. At least one carbon atom that is carbonyl (double-bonded to an O).
    for ring in ring_info:
        if len(ring) == 6:
            # Check if any atom in the ring is oxygen.
            has_ring_oxygen = any(mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring)
            
            # Check if any carbon in the ring has a double bond to an oxygen (a carbonyl group).
            has_carbonyl = False
            for i in ring:
                atom = mol.GetAtomWithIdx(i)
                # Consider carbons only
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Check bond type (double bond). (RDKit bond types compared to 2 are double bonds.)
                        if bond.GetBondTypeAsDouble() == 2:
                            # See if the other atom is oxygen
                            other = bond.GetOtherAtom(atom)
                            if other.GetAtomicNum() == 8:
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            
            if has_ring_oxygen and has_carbonyl:
                candidate_rings.append(set(ring))
    
    # Look for a pair of candidate rings that are fused (i.e. share at least 2 atoms).
    fused_pair_found = False
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            # If the intersection contains two or more atoms, we consider these rings fused.
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                fused_pair_found = True
                break
        if fused_pair_found:
            break

    if fused_pair_found:
        return True, ("Found fused six‐membered ring system with a chromenone-like ring fused to an "
                      "oxygen-containing ring, consistent with a rotenoid skeleton")
    else:
        return False, "No characteristic fused chromenone ring system found; molecule is unlikely to be a rotenoid"

# Example test(s) – you may run these (although in production a full suite would be used)
if __name__ == "__main__":
    # Example: 13alpha-Hydroxydolineone (one of the given rotenoids)
    test_smiles = "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C=C6OC=CC6=C5"
    result, reason = is_rotenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)