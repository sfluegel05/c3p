"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty acyl-CoA(4-)
Definition: A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
of any long-chain fatty acyl-CoA; major species at pH 7.3.
This implementation verifies that the molecule contains (a) an adenine nucleotide fragment (using three alternative SMARTS patterns to allow for substituents),
(b) a pantetheine fragment, (c) a thioester group connecting a fatty acyl chain and (d) a linear carbon chain (12–22 carbons) attached to the thioester.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) belongs to the class
    long-chain fatty acyl-CoA(4-).

    Steps:
      1. Parse the molecule.
      2. Verify it is a CoA derivative by checking for:
         a. An adenine nucleotide fragment. Three SMARTS patterns are attempted to catch varying representations.
         b. The pantetheine fragment.
      3. Look for a thioester group (a non-ring carbonyl carbon double-bonded to O and bonded to S).
      4. For each found thioester, follow from the carbonyl carbon (excluding the S neighbor) into a linear (unbranched) chain.
      5. Confirm that one of these chains is 12–22 carbons long.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as long-chain fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the adenine nucleotide fragment using three alternative SMARTS patterns.
    adenine_pattern1 = Chem.MolFromSmarts("n1cnc2ncnc12")           # no substituents
    adenine_pattern2 = Chem.MolFromSmarts("n1cnc2c(n)ncnc12")         # allowing lower-case substituent on one ring atom
    adenine_pattern3 = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")         # explicit NH2 substituent variant
    has_adenine = (mol.HasSubstructMatch(adenine_pattern1) or 
                   mol.HasSubstructMatch(adenine_pattern2) or 
                   mol.HasSubstructMatch(adenine_pattern3))
    if not has_adenine:
        return False, "Adenine nucleotide fragment not found; not a CoA derivative"
    
    # Check for the pantetheine fragment characteristic of CoA.
    pant_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pant_pattern):
        return False, "Pantetheine fragment not found; not a CoA derivative"
    
    # Look for a thioester group. Pattern: a non-ring carbon (C;!R) double-bonded to O and single-bonded to S.
    thioester_pattern = Chem.MolFromSmarts("[C;!R](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (fatty acid linkage) found"
    
    # Helper function: recurisvely count contiguous, unbranched carbon atoms.
    def linear_chain_length(atom, from_idx):
        """
        Recursively counts contiguous carbon atoms (unbranched) from the given atom.
        The recursion ends if there is not exactly one carbon neighbor apart from the atom where we came.
        """
        length = 1  # count the current atom
        # Filter carbon neighbors, excluding the one we came from (from_idx)
        nbr_carbons = [nbr for nbr in atom.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != from_idx]
        # If chain branches (or dead ends), return the current chain count.
        if len(nbr_carbons) != 1:
            return length
        next_atom = nbr_carbons[0]
        # Check the next step: if the next atom (beyond the one we came from) has >1 carbon neighbor,
        # count it and stop.
        further_nbrs = [nbr for nbr in next_atom.GetNeighbors()
                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
        if len(further_nbrs) > 1:
            return length + 1
        return length + linear_chain_length(next_atom, atom.GetIdx())
    
    # Specify acceptable carbon count limits for the fatty chain.
    MIN_CARBONS = 12
    MAX_CARBONS = 22
    acyl_chain_found = False
    found_chain_length = 0
    
    # For every thioester match, try to trace the fatty acyl chain.
    for match in thioester_matches:
        carbonyl_idx = match[0]  # carbonyl carbon index
        sulfur_idx = match[1]     # sulfur index
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From the carbonyl atom, choose a neighbor that is not the sulfur and is carbon
        fatty_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_start = nbr
                break
        if fatty_start is None:
            continue  # try next thioester match if no candidate fatty acid carbon is found
        # Count the linear chain length from the fatty acid starting atom.
        chain_length = linear_chain_length(fatty_start, carbonyl_idx)
        if MIN_CARBONS <= chain_length <= MAX_CARBONS:
            acyl_chain_found = True
            found_chain_length = chain_length
            break
    
    if not acyl_chain_found:
        return False, (f"Fatty acyl chain missing or not linear with required {MIN_CARBONS}-{MAX_CARBONS} carbons.")
    
    # All tests passed: the molecule contains an adenine nucleotide, pantetheine fragment,
    # a thioester group, and an acceptable fatty acyl chain.
    return True, (
        f"Contains CoA moiety (adenine and pantetheine fragments), a thioester group with a fatty acyl chain of length {found_chain_length} carbons, "
        "and expected connectivity for long-chain fatty acyl-CoA(4-)."
    )

# Example usage (for testing):
if __name__ == "__main__":
    # Test the example with one of the provided SMILES strings.
    test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)