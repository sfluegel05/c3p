"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar - any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
"""

from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is defined as any monosaccharide that contains an alcoholic hydroxy group 
    that is esterified with phosphoric acid (i.e. an -O-P(=O)(O)O moiety attached to a sugar ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a phospho sugar, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify candidate sugar rings.
    # We consider a ring of 5 atoms (furanose) or 6 atoms (pyranose)
    # that contains exactly one oxygen atom (the ring-oxygen) as a typical monosaccharide.
    ri = mol.GetRingInfo()
    candidate_sugar_atom_indices = set()  # indices of carbon atoms that are part of a candidate sugar ring
    has_candidate_ring = False
    for ring in ri.AtomRings():
        if len(ring) not in {5, 6}:
            continue
        # Count number of oxygen atoms in the ring.
        oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
        # For a typical sugar ring, we expect exactly one ring oxygen.
        if oxygen_count == 1:
            has_candidate_ring = True
            # Add the carbon atoms in the ring as candidate points for the phosphate ester.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon
                    candidate_sugar_atom_indices.add(idx)
    
    if not has_candidate_ring:
        return False, "No appropriate monosaccharide ring (5- or 6-membered with one oxygen) found"
    
    # Now, look for an alcoholic -OH attached to one of the candidate sugar carbons that is esterified with phosphoric acid.
    # We require: sugar carbon - O - P, with P having at least one double-bonded oxygen.
    for cidx in candidate_sugar_atom_indices:
        carbon = mol.GetAtomWithIdx(cidx)
        for nbr in carbon.GetNeighbors():
            # Look for an oxygen neighbor (exocyclic hydroxyl group)
            if nbr.GetAtomicNum() != 8:
                continue
            # If this oxygen is part of the ring, skip it (we want the alcoholic group that is exocyclic)
            if nbr.GetIdx() in candidate_sugar_atom_indices:
                continue
            # For this oxygen, check its neighbors (aside from our sugar carbon).
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == cidx:
                    continue
                # Look for a phosphorus (atomic num 15)
                if nbr2.GetAtomicNum() == 15:
                    phosphorus = nbr2
                    # Check that phosphorus has at least one double-bonded oxygen.
                    has_dbl_bonded_oxygen = False
                    for p_nbr in phosphorus.GetNeighbors():
                        # Skip the oxygen that connected back to sugar.
                        if p_nbr.GetAtomicNum() != 8:
                            continue
                        # Retrieve the bond between phosphorus and this oxygen.
                        bond = mol.GetBondBetweenAtoms(phosphorus.GetIdx(), p_nbr.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                            has_dbl_bonded_oxygen = True
                            break
                    if has_dbl_bonded_oxygen:
                        return True, "Found sugar ring with an alcoholic –OH esterified to a phosphoric acid"
    
    return False, "No alcoholic –OH on a candidate sugar ring found that is esterified to phosphoric acid"

# For testing, you could run:
if __name__ == "__main__":
    # Example: alpha-D-Hexose 6-phosphate
    test_smiles = "P(OCC1O[C@H](O)C(O)C(O)C1O)(O)(O)=O"
    result, reason = is_phospho_sugar(test_smiles)
    print(f"Is phospho sugar? {result}\nReason: {reason}")