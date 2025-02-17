"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose (in our definition) is a hexose (exactly 6 carbons and 6 oxygens) that contains a sugar ring 
(5- or 6-membered with exactly one oxygen) and has exactly one ring-carbon that carries an exocyclic –CH2OH 
group (the “C5” center). That C5 center must have CIP configuration 'R'.
"""

from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines whether a molecule is a D-hexose (D-configured aldose) based on its SMILES string.
    
    A D-hexose here is required to:
      - Have exactly 6 carbon atoms and 6 oxygen atoms.
      - Contain a sugar ring that is 5- or 6-membered and contains exactly one oxygen.
      - Have exactly one exocyclic CH2OH substituent attached to a ring carbon. This ring carbon (the C5
        center) must have CIP configuration 'R'.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        (bool, str): Tuple of classification result and an explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that CH2OH groups have their hydrogen count explicit.
    mol = Chem.AddHs(mol)
    
    # Assign stereochemistry so that CIP codes are computed.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Count heavy atoms: exactly 6 carbons and exactly 6 oxygens are required.
    atoms = list(mol.GetAtoms())
    carbons = [a for a in atoms if a.GetAtomicNum() == 6]
    oxygens = [a for a in atoms if a.GetAtomicNum() == 8]
    if len(carbons) != 6:
        return False, "Molecule does not have exactly 6 carbon atoms (not a hexose)"
    if len(oxygens) != 6:
        return False, "Molecule does not have exactly 6 oxygen atoms (likely a modified hexose)"
    
    # Get ring information and look for a sugar ring: a 5- or 6-membered ring that contains exactly one oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_ring = None
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms in the ring.
        ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if ring_oxygens != 1:
            continue
        # We found a candidate sugar ring.
        candidate_ring = set(ring)
        break
    if candidate_ring is None:
        return False, "No appropriate sugar ring (5- or 6-membered with one oxygen) found"
    
    # Now find a candidate exocyclic CH2OH group attached to one of the ring carbons.
    # We define a valid CH2OH candidate as:
    #   - an sp3 carbon (atomic number 6) that is NOT in the ring,
    #   - has exactly 2 hydrogens attached (GetTotalNumHs() == 2),
    #   - is attached (via a single bond) to exactly one oxygen (which should be -OH).
    exo_candidates = []
    # Loop over atoms in the ring. Only consider carbons.
    for idx in candidate_ring:
        ring_atom = mol.GetAtomWithIdx(idx)
        if ring_atom.GetAtomicNum() != 6:
            continue
        # Check among neighbors which are NOT in the candidate ring.
        for nbr in ring_atom.GetNeighbors():
            if nbr.GetIdx() in candidate_ring:
                continue
            # We are looking for a CH2OH substituent: first the neighbor must be carbon.
            if nbr.GetAtomicNum() != 6:
                continue
            # Check for exactly 2 attached hydrogens (explicit + implicit)
            if nbr.GetTotalNumHs() != 2:
                continue
            # Now, candidate nbr should be attached to exactly one oxygen (forming the -OH) among its neighbors,
            # aside from the bond back to the ring.
            oxy_neighbors = []
            for subnbr in nbr.GetNeighbors():
                # Exclude the ring_atom we already came from.
                if subnbr.GetIdx() == ring_atom.GetIdx():
                    continue
                if subnbr.GetAtomicNum() == 8:
                    # Make sure the oxygen is attached via a single bond.
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 1:
                        oxy_neighbors.append(subnbr)
            if len(oxy_neighbors) != 1:
                continue
            # We now have a candidate CH2OH substituent.
            exo_candidates.append((ring_atom, nbr))
    
    if len(exo_candidates) == 0:
        return False, "No appropriate exocyclic CH2OH substituent found"
    if len(exo_candidates) > 1:
        return False, "Multiple exocyclic CH2OH substituents found (ambiguous candidate for C5)"
    
    # Exactly one candidate was found. The ring atom to which the exocyclic CH2OH is attached is our candidate C5.
    candidate_C5, exo_substituent = exo_candidates[0]
    
    # Check that the candidate_C5 has valid stereochemistry assigned and that it is 'R'
    if not candidate_C5.HasProp('_CIPCode'):
        return False, "Stereochemistry was not assigned for the candidate C5 atom"
    cip = candidate_C5.GetProp('_CIPCode')
    if cip != 'R':
        return False, f"C5 has CIP configuration {cip} (expected R for D-hexose)"
    
    return True, "Contains a sugar ring (5- or 6-membered with one oxygen) with exactly one exocyclic CH2OH substituent attached to a ring-carbon (C5) having R configuration"

# Example usage and basic testing:
if __name__ == "__main__":
    test_smiles_list = [
        # True positives
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)CO",             # beta-D-idopyranose
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",                   # alpha-D-galactose
        # False positives (should be rejected)
        "OC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",                     # beta-D-fructopyranose (ketose/lactone-like)
        "OC[C@H]1OC(=O)[C@H](O)[C@@H](O)[C@H]1O",                       # D-galactono-1,5-lactone
        # False negatives (should be accepted but were missed before but may still have issues)
        "O1[C@H]([C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)CO"                # beta-D-idofuranose
    ]
    for smi in test_smiles_list:
        res, reason = is_D_hexose(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*40}")