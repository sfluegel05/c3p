"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
#!/usr/bin/env python3
"""
Classifies: Tertiary amine oxide
Definition: An N-oxide in which a tetrahedral (sp³) nitrogen (with a +1 formal charge) 
is bonded to exactly one oxygen atom (with a -1 formal charge via a single bond) and three carbon atoms.
Furthermore, for each carbon substituent the atoms immediately bonded (except the candidate N) must be only carbons or hydrogens,
unless the bond to that neighbor is part of a ring that is composed exclusively of carbons and hydrogens.
"""

from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule contains a tertiary amine oxide functional group based on its SMILES string.
    
    Improved algorithm:
      1. Parse the SMILES string.
      2. Look for candidate nitrogen atoms with:
         - Symbol "N" and formal charge of +1.
         - Tetrahedral geometry (degree = 4) and not aromatic.
      3. For each candidate N, require exactly one neighbor oxygen atom (with -1 formal charge) 
         that is single-bonded and has no other attachment (degree == 1) and exactly three carbon neighbors.
      4. For each carbon neighbor, check every immediate neighbor (except the candidate N). 
         For each such bond:
           - If the neighbor atom is carbon or hydrogen, it is allowed.
           - Otherwise, if the bond between the carbon and that atom is part of a ring AND at least one ring 
             that contains both atoms is composed only of carbons and hydrogens, then the connection is allowed.
           - Otherwise, the substituent is too complex.
      5. If a candidate nitrogen satisfies all the tests, return True and a clear explanation.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a tertiary amine oxide functional group is found, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information: each ring is reported as a tuple of atom indices.
    rings = mol.GetRingInfo().AtomRings()

    def bond_in_pure_CH_ring(idx1, idx2):
        """
        Checks if the bond between atom indices idx1 and idx2 is part of at least one ring that is composed solely of C and H.
        """
        # Obtain the bond between idx1 and idx2.
        bond = mol.GetBondBetweenAtoms(idx1, idx2)
        if bond is None:
            return False
        # Only consider bonds that are in a ring.
        if not bond.IsInRing():
            return False
        # Iterate through rings that contain both atoms.
        for ring in rings:
            if idx1 in ring and idx2 in ring:
                # Check if every atom in this ring is either carbon (atomic num 6) or hydrogen (atomic num 1)
                if all(mol.GetAtomWithIdx(a).GetAtomicNum() in (1,6) for a in ring):
                    return True
        return False

    # For each potential candidate nitrogen.
    for atom in mol.GetAtoms():
        # Check symbol and formal charge.
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 1:
            continue
        # Require tetrahedral geometry (degree==4) and not aromatic.
        if atom.GetDegree() != 4 or atom.GetIsAromatic():
            continue
        
        neighbors = atom.GetNeighbors()
        oxide_count = 0
        candidate_oxygen = None
        carbon_neighbors = []
        
        # Examine the neighbors.
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Look for the oxygen substituent with -1 charge and a single bond.
            if nbr.GetSymbol() == "O" and nbr.GetFormalCharge() == -1:
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxide_count += 1
                    candidate_oxygen = nbr
            elif nbr.GetAtomicNum() == 6:  # Must be a carbon atom.
                carbon_neighbors.append(nbr)
            else:
                # Any other substituents are not acceptable.
                pass
        
        # Must have exactly one oxide and three carbon substituents.
        if oxide_count != 1 or len(carbon_neighbors) != 3:
            continue
        
        # Check that the oxide oxygen is bonded only to this candidate N.
        if candidate_oxygen.GetDegree() != 1:
            continue

        # For each carbon neighbor, check its immediate neighbors (except back to the candidate N).
        substituents_ok = True
        for c_atom in carbon_neighbors:
            for nbr in c_atom.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue  # Skip the bond back to the candidate nitrogen.
                # If neighbor is carbon or hydrogen, allowed.
                if nbr.GetAtomicNum() in (1,6):
                    continue
                # Otherwise, require that the bond from c_atom to nbr is part of a pure CH ring.
                if not bond_in_pure_CH_ring(c_atom.GetIdx(), nbr.GetIdx()):
                    substituents_ok = False
                    break
            if not substituents_ok:
                break

        if not substituents_ok:
            continue  # This candidate N fails the substituent simplicity test.
        
        # Found a valid tertiary amine oxide functional group.
        return True, ("Found a tertiary amine oxide functional group "
                      "(N+ bonded to O- via a single bond and three carbon substituents, each of which is a simple organic group).")
    
    # No candidate passed all tests.
    return False, "No tertiary amine oxide functional group found"


# Optional main testing block.
if __name__ == "__main__":
    # Some example SMILES from the provided outcomes (both expected positives and negatives)
    test_examples = [
        # True positives:
        ("[H][C@]1(C[N+]2([O-])CC[C@@]34C(Nc5ccccc35)=C(C(=O)OC)[C@@]1([H])C[C@]24[H])C(C)=O", "(-)-Alstolucine C"),
        ("C[N+](C)([O-])C", "trimethylamine N-oxide"),
        ("C1=CC([C@]2([N@+](CCC2)([O-])C)[H])=CN=C1", "trans-(S)-nicotine N(1')-oxide"),
        ("C[N+]1([O-])CCC(=CC1)c1ccccc1", "1-methyl-4-phenyl-1,2,3,6-tetrahydropyridine N-oxide"),
        # False positive examples (should report negative):
        ("C(CC[N+](C)(C)[O-])N1C=2C(CCC3=C1C=CC=C3)=CC=CC2", "imipramine oxide (should not be classified)"),
        ("O[C@]1([C@H](CCCC1)C[N+]([O-])(C)C)C2=CC(OC)=CC=C2", "Tramadol N-Oxide (should not be classified)"),
    ]
    
    for s, name in test_examples:
        result, explanation = is_tertiary_amine_oxide(s)
        print(f"SMILES: {s}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")