"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
#!/usr/bin/env python
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
The molecule must contain a carnitine moiety – a chiral carbon attached to 
a carboxylate group (CC([O-])=O) and a trimethylammonium substituent (C[N+](C)(C)C)
– and its hydroxyl group must be acylated (i.e. bound to a carbonyl via an ester linkage).
We first use a non-chiral SMARTS to find a carnitine-like core; then we assign CIP labels 
and check that the carnitine stereocenter is (R)-configured (which is equivalent to L-carnitine). 
Finally, we verify that the oxygen in the carnitine core is acylated by a carbonyl group.
"""

from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.

    The procedure is:
      1. Parse the SMILES and assign stereochemistry.
      2. Look for a carnitine-like substructure (without requiring chirality markers) 
         using the pattern "O[C](CC([O-])=O)C[N+](C)(C)C". In this pattern:
             - The oxygen (first atom) is the site where the acyl ester forms.
             - The second atom (C) is the carnitine chiral center.
             - The third atom (from "CC([O-])=O") is the carboxylate side.
             - The fourth fragment is the trimethylammonium group.
      3. For each matching instance, retrieve the chiral center from the carnitine core
         and use RDKit’s stereochemistry tools to get its CIP label – for L-carnitine we expect "R".
      4. Then verify that the oxygen in the core (the first atom in the match) is bound
         to an acyl group, i.e. one of its neighbors (other than the carnitine center) is a carbon
         that is part of a carbonyl group (C=O).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as an O-acyl-L-carnitine, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the input SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned (this computes CIP labels)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define a non-chiral SMARTS pattern to find the carnitine-like core.
    # Pattern: an oxygen connected to a carbon (the chiral center) that is also bound to:
    #   - a carboxylate fragment: CC([O-])=O, and
    #   - a trimethylammonium group: C[N+](C)(C)C.
    # We do not impose chirality in the SMARTS so that we can verify it later via CIP.
    carnitine_pattern = Chem.MolFromSmarts("O[C](CC([O-])=O)C[N+](C)(C)C")
    matches = mol.GetSubstructMatches(carnitine_pattern)
    if not matches:
        return False, "No carnitine-like substructure found."

    # Process each match.
    for match in matches:
        # In our pattern, we expect:
        #    match[0] = the oxygen attached to the acyl group,
        #    match[1] = the carnitine chiral carbon,
        #    match[2] = part of the carboxylate group,
        #    match[3] = the N (or part of the trimethylammonium group).
        chiral_idx = match[1]
        chiral_atom = mol.GetAtomWithIdx(chiral_idx)
        # Verify that a CIP label (R/S) has been assigned.
        try:
            cip = chiral_atom.GetProp('_CIPCode')
        except KeyError:
            cip = None
        # For L-carnitine, the chiral center (in its proper depiction) must be (R).
        if cip != "R":
            # Skip this match if the configuration is not "R".
            continue
        
        # Check that the oxygen atom (match[0]) is acylated.
        o_atom = mol.GetAtomWithIdx(match[0])
        neighbors = o_atom.GetNeighbors()
        # We expect 2 neighbors: one is the carnitine chiral carbon, and the other should be an acyl carbon.
        # Remove the carnitine carbon.
        acyl_neighbor = None
        for nb in neighbors:
            if nb.GetIdx() != chiral_idx:
                acyl_neighbor = nb
                break
        if acyl_neighbor is None:
            continue
        
        # The acyl neighbor should be a carbon.
        if acyl_neighbor.GetAtomicNum() != 6:
            continue
        
        # Check if the acyl neighbor is part of a carbonyl group (i.e. it has at least one double-bonded O neighbor).
        carbonyl_found = False
        for nb in acyl_neighbor.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(acyl_neighbor.GetIdx(), nb.GetIdx())
            if nb.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                carbonyl_found = True
                break
        if not carbonyl_found:
            continue
        
        # If all checks pass, we have a carnitine core with the correct L (i.e. R CIP) configuration and proper acylation.
        return True, "Molecule contains an O-acyl-L-carnitine motif with proper acylation and L-carnitine configuration."
    
    return False, "Carnitine-like substructure found, but either configuration is not L (R CIP) or the acyl group attachment is missing."
    
# Example usage:
if __name__ == "__main__":
    # A list of example SMILES strings, including some known positives.
    test_smiles = [
        "CCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-butanoyl-L-carnitine: should be True.
        "CC(=O)O[C@@H](CC([O-])=O)C[N+](C)(C)C"      # O-acetyl-D-carnitine: configuration is not (R) so should be False.
    ]
    for s in test_smiles:
        result, reason = is_O_acyl_L_carnitine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")