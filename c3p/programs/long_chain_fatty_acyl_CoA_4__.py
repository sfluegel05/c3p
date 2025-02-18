"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty acyl-CoA(4-)
Definition: A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
of any long-chain fatty acyl-CoA; major species at pH 7.3.
This version uses two alternative adenine SMARTS patterns (to cover possible representations),
the pantetheine fragment, and a thioester group search along with a strict linear-chain search
to identify a fatty acyl chain of 12–22 contiguous carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) belongs to the class
    long-chain fatty acyl-CoA(4-).

    The algorithm is as follows:
      1. Parse the molecule.
      2. To verify it is a CoA derivative, check for:
         a. A nucleotide (adenine) fragment. Two SMARTS patterns are attempted.
         b. The pantetheine fragment.
      3. Locate a thioester group (a non-ring carbonyl carbon double-bonded to oxygen and single-bonded to a sulfur).
      4. For each thioester match, follow the carbon attached (other than the sulfur) using a helper
         function that “walks” a strictly linear (i.e. unbranched) carbon chain.
      5. Accept the molecule if a chain between 12 and 22 carbons is found.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as long-chain fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Instead of using an overall charge, we verify that two key fragments of CoA are present.

    # Check for the adenine nucleotide fragment by trying two alternative adenine SMARTS patterns.
    adenine_pattern1 = Chem.MolFromSmarts("n1cnc2ncnc12")
    adenine_pattern2 = Chem.MolFromSmarts("n1cnc2c(n)ncnc12")
    has_adenine = mol.HasSubstructMatch(adenine_pattern1) or mol.HasSubstructMatch(adenine_pattern2)
    if not has_adenine:
        return False, "Adenine nucleotide fragment not found; not a CoA derivative"

    # Check for the pantetheine fragment characteristic of CoA.
    pant_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pant_pattern):
        return False, "Pantetheine fragment not found; not a CoA derivative"

    # Look for a thioester group. The pattern looks for a non-ring carbonyl carbon double bonded to O and single bonded to S.
    thioester_pattern = Chem.MolFromSmarts("[C;!R](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (fatty acid linkage) found"

    # Define a helper function that follows a strictly linear, unbranched carbon chain.
    def linear_chain_length(atom, from_idx):
        """
        Recursively counts contiguous carbon atoms from 'atom' excluding coming back to from_idx.
        If the chain branches (has more than one carbon neighbor besides the one already visited),
        the walk stops.
        """
        length = 1  # count the current atom
        # Get carbon neighbors excluding the one we came from
        nbr_carbons = [nbr for nbr in atom.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != from_idx]
        # If branching or no continuation, return current length.
        if len(nbr_carbons) != 1:
            return length
        next_atom = nbr_carbons[0]
        # Check linearity: if next_atom (beyond the coming atom) has more than one carbon neighbor, stop after counting it.
        further_nbrs = [nbr for nbr in next_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
        if len(further_nbrs) > 1:
            return length + 1
        return length + linear_chain_length(next_atom, atom.GetIdx())

    # Define the acceptable length boundaries for a long-chain fatty acyl chain.
    MIN_CARBONS = 12
    MAX_CARBONS = 22

    acyl_chain_found = False
    found_chain_length = 0

    # The thioester pattern returns two atoms (carbonyl carbon and sulfur).
    for match in thioester_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From the carbonyl atom, find a neighboring carbon (excluding the sulfur) that is likely 
        # the beginning of the fatty acyl chain.
        fatty_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_start = nbr
                break
        if fatty_start is None:
            continue  # try the next thioester match if no candidate fatty carbon is found
        # Count the linear carbon chain length starting from fatty_start.
        chain_length = linear_chain_length(fatty_start, carbonyl_idx)
        if MIN_CARBONS <= chain_length <= MAX_CARBONS:
            acyl_chain_found = True
            found_chain_length = chain_length
            break

    if not acyl_chain_found:
        return False, (f"Fatty acyl chain missing or not linear and within {MIN_CARBONS}-{MAX_CARBONS} carbons.")

    # Classification successful: CoA moiety (adenine and pantetheine fragments) and an acceptable thioester-linked fatty chain.
    return True, (
        f"Contains CoA moiety (adenine and pantetheine fragments), "
        f"a thioester group with a fatty acyl chain of length {found_chain_length} carbons, "
        "and expected connectivity for long-chain fatty acyl-CoA(4-)."
    )

# Example usage:
if __name__ == "__main__":
    # Test with one of the positive examples.
    test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)