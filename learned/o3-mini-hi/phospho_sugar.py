"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar - any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
The definition used here is that a molecule must contain a sugar ring (a furanose or pyranose substructure)
and an exocyclic oxygen (an –OH group) on that ring that is connected (via a single bond) to a phosphorus atom,
and that phosphorus must carry at least one double‐bonded oxygen.
"""

from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is defined as any monosaccharide that contains an alcoholic -OH group 
    that is esterified with phosphoric acid (i.e. an -O-P(=O)(O)O moiety attached to a sugar ring).
    
    This implementation:
    1. Attempts to identify a sugar ring by matching typical pyranose or furanose SMARTS patterns.
    2. For each sugar ring found, it scans the ring atoms and collects any exocyclic oxygen.
    3. If any such oxygen is directly bound (by a single bond) to a phosphorus atom whose bonding pattern
       includes at least one double-bonded oxygen, the function returns True.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a phospho sugar, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for typical sugar rings.
    # Pyranose: six membered ring with five carbons and one oxygen.
    pyranose_smarts = "[C;R]1[C;R][C;R][C;R][C;R][O;R]1"
    # Furanose: five membered ring with four carbons and one oxygen.
    furanose_smarts = "[C;R]1[C;R][C;R][C;R][O;R]1"
    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    furanose = Chem.MolFromSmarts(furanose_smarts)
    
    sugar_ring_matches = []
    
    if mol.HasSubstructMatch(pyranose):
        sugar_ring_matches.extend(mol.GetSubstructMatches(pyranose))
    if mol.HasSubstructMatch(furanose):
        sugar_ring_matches.extend(mol.GetSubstructMatches(furanose))
    
    if not sugar_ring_matches:
        return False, "No sugar ring (pyranose or furanose pattern) found"
    
    # Helper function to check if an oxygen (attached to the sugar ring) is part of a phosphate ester.
    def oxygen_has_phosphate_connection(oxygen, mol):
        # Look at neighbors of the oxygen.
        for nb in oxygen.GetNeighbors():
            if nb.GetAtomicNum() == 15:  # phosphorus
                # Check that the O-P bond is single.
                bond_op = mol.GetBondBetweenAtoms(oxygen.GetIdx(), nb.GetIdx())
                if bond_op is None or bond_op.GetBondType() != Chem.BondType.SINGLE:
                    continue
                phosphorus = nb
                has_dbl_bonded_O = False
                # Check phosphorus neighbors: it should have at least one O in a double bond.
                for p_nb in phosphorus.GetNeighbors():
                    if p_nb.GetAtomicNum() != 8:  # oxygen
                        continue
                    bond_po = mol.GetBondBetweenAtoms(phosphorus.GetIdx(), p_nb.GetIdx())
                    if bond_po is not None and bond_po.GetBondType() == Chem.BondType.DOUBLE:
                        has_dbl_bonded_O = True
                        break
                if has_dbl_bonded_O:
                    return True
        return False

    # For each sugar ring match, look for an exocyclic oxygen connected to the ring.
    for ring in sugar_ring_matches:
        ring_atom_indices = set(ring)  # indices that belong to the sugar ring
        # Flag for found phosphate ester on this sugar ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # For each neighbor of the ring atom:
            for nbr in atom.GetNeighbors():
                # Skip if the neighbor is also in the ring (we want exocyclic substituents)
                if nbr.GetIdx() in ring_atom_indices:
                    continue
                # We are interested in oxygen atoms (potential -OH group)
                if nbr.GetAtomicNum() != 8:
                    continue
                # Optionally, one could check if the oxygen is sp3-hybridized.
                if nbr.GetHybridization() != Chem.HybridizationType.SP3:
                    continue
                # Now, check if this oxygen is bound to a phosphate group (P with a =O)
                if oxygen_has_phosphate_connection(nbr, mol):
                    return True, "Found sugar ring with an alcoholic –OH esterified to a phosphoric acid"
    
    return False, "No alcoholic –OH on a candidate sugar ring found that is esterified with phosphoric acid"

# For testing (run if executed as a script):
if __name__ == "__main__":
    # Example test: alpha-D-Hexose 6-phosphate
    test_smiles = "P(OCC1O[C@H](O)C(O)C(O)C1O)(O)(O)=O"
    result, reason = is_phospho_sugar(test_smiles)
    print(f"SMILES: {test_smiles}\nIs phospho sugar? {result}\nReason: {reason}")