"""
Classifies: CHEBI:17478 aldehyde
"""
#!/usr/bin/env python
"""
Classifies: Aldehyde – a compound having the functional group RC(=O)H 
(including formaldehyde H2C=O) where the carbonyl carbon is bonded 
to exactly one hydrogen and one R group (or two hydrogens in formaldehyde).
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a compound containing a carbonyl group (C=O)
    in which the carbonyl carbon is bonded to exactly one hydrogen and one R group,
    or in the case of formaldehyde, two hydrogens.
    
    We first add explicit hydrogens to ensure connectivity is reliable,
    then search for carbonyl carbons with the SMARTS "[CX3](=O)".
    For each candidate, we:
     1. Check that the carbon has exactly three neighbors.
     2. Identify the double-bonded oxygen (by checking for a bond order of 2).
     3. Check that the remaining two neighbors are either both hydrogens
        (formaldehyde) or one hydrogen and one non-oxygen (typical aldehyde).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if a valid aldehyde functional group is found, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to reliably inspect connectivity.
    mol = Chem.AddHs(mol)

    # Define a SMARTS pattern to match any carbonyl carbon.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)")
    if carbonyl_pattern is None:
        return False, "Failed to create carbonyl pattern"
    
    # Get all substructure matches for carbonyl carbons.
    candidate_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not candidate_matches:
        return False, "No carbonyl group (C=O) found in the molecule"
    
    valid_count = 0
    # Iterate through each candidate carbonyl carbon.
    for match in candidate_matches:
        carbon = mol.GetAtomWithIdx(match[0])
        # Get all neighbors of the carbon atom.
        neighbors = carbon.GetNeighbors()
        # For an aldehyde, the carbon (after adding explicit Hs) should be bound to 3 atoms:
        # one is the double-bonded oxygen; the other two are substituents.
        if len(neighbors) != 3:
            continue

        dO = None  # to store the double-bonded oxygen atom
        other_neighbors = []  # list for the remaining two substituents
        
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond is None:
                continue
            # Check if the neighbor is double-bonded oxygen.
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                # There should be exactly one double-bonded oxygen.
                if dO is not None:
                    dO = None
                    break
                dO = nbr
            else:
                other_neighbors.append(nbr)
        
        # Ensure exactly one carbonyl oxygen is found.
        if dO is None:
            continue

        # Check that there are exactly two other substituents.
        if len(other_neighbors) != 2:
            continue

        # Evaluate the two remaining substituents:
        # Count how many are hydrogens and record the non-hydrogen (if any).
        h_count = 0
        non_h_neighbor = None
        for nbr in other_neighbors:
            if nbr.GetAtomicNum() == 1:
                h_count += 1
            else:
                non_h_neighbor = nbr

        # Valid connectivity conditions:
        # 1) Formaldehyde: both neighbors are hydrogens (h_count == 2).
        # 2) Typical aldehyde: one hydrogen and one non-hydrogen (h_count == 1)
        #    and the non-hydrogen must not be oxygen.
        if (h_count == 2) or (h_count == 1 and non_h_neighbor is not None and non_h_neighbor.GetAtomicNum() != 8):
            valid_count += 1

    if valid_count == 0:
        return False, "No valid aldehyde group found after checking connectivity"
    
    reason = f"Found {valid_count} aldehyde group{'s' if valid_count > 1 else ''} with proper connectivity (C=O with one carbonyl O, one H and one R group or H2 for formaldehyde)"
    return True, reason

# Optional testing code; remove or comment out when using in a larger project.
if __name__ == "__main__":
    test_smiles = [
        "O=CC(CCC=C(C)C)C",  # 5-Heptenal, 2,6-dimethyl-
        "Oc1c(C=O)ccc2ccccc12",  # 1-hydroxy-2-naphthaldehyde
        "CC(=O)O[C@H]1CC[C@]2(C=O)[C@H]3CC[C@@]4(C)[C@@H](CCC4=O)[C@@H]3CC=C2C1",  # 3beta-Hydroxy-17-oxoandrost-5-en-19-al acetate
        "C(CCCCCCCCCCC(=O)[H])CCCCC",  # heptadecanal
        "O=CCCCCCCCCCCCCCCCCCCCCCCC",  # tetracosanal
        "CCCC/C=C/C=O",  # (E)-hept-2-enal
        "CCCCCCCCCCCCCCCCC=O",  # octadecanal
        "[H]C(=O)c1cc(N)ccc1O",  # 5-aminosalicylaldehyde
        "[H]C(=O)c1cccc(O)c1O",  # 2,3-dihydroxybenzaldehyde
        "Cc1cc(O)c(C=O)c2Oc3cc4OC(C)(C)CC(=O)c4c(C)c3OC(=O)c12",  # mollicellin B
        "C[C@H](CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C)C=O",  # (25R)-3-ketocholest-4-en-26-al
        "OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC=O",  # 20-oxoarachidonic acid
        "C(C#C)=O",  # prop-2-ynal should classify as aldehyde
        "[H]C(=O)",  # formaldehyde with explicit hydrogen
        "CC(C)CC=O",  # 3-methylbutanal
    ]
    
    for sm in test_smiles:
        result, explanation = is_aldehyde(sm)
        print(f"SMILES: {sm}\nResult: {result}\nExplanation: {explanation}\n")