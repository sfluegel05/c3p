"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty acyl-CoA(4-)
Definition: A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
of any long-chain fatty acyl-CoA; major species at pH 7.3.
This updated program uses both substructure searches for a nucleotide fragment and the pantetheine part
to reliably detect a CoA moiety. It also looks for a thioester group linking a fatty acyl chain and then
counts the number of contiguous (and unbranched) carbon atoms. Only molecules with a linear chain of 
between 12 and 22 carbons (inclusive) are accepted as long‐chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) belongs to the class
    long-chain fatty acyl-CoA(4-).

    The algorithm now does the following:
      1. Parses the molecule.
      2. Rather than relying on overall formal charge, it checks that the molecule contains
         two key pieces of a CoA: (a) a nucleotide fragment (adenine part) and (b) the pantetheine fragment.
      3. It locates a thioester group (pattern: a carbonyl carbon bound to a sulfur).
      4. For each thioester match it identifies a candidate fatty acyl chain by following the carbon
         attached (other than the sulfur) to the carbonyl. It then walks only through carbon atoms 
         in a strictly linear, unbranched fashion. (Any deviation – if a carbon has more than one carbon neighbor 
         that isn’t the one we came from – will abort that path.)
      5. The chain length must be between 12 and 22 contiguous carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as long-chain fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Instead of checking the overall charge (which may be omitted in SMILES), we ensure that the CoA unit is present.
    # Check for the adenine nucleotide fragment common in CoA.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine nucleotide fragment not found; not a CoA derivative"

    # Check for the pantetheine fragment part – this fragment is characteristic of CoA.
    pant_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pant_pattern):
        return False, "Pantetheine fragment not found; not a CoA derivative"

    # Look for a thioester group: a carbonyl carbon (not in a ring) double bonded to O and single bonded to S.
    # The SMARTS "[C;!R](=O)[S]" is used.
    thioester_pattern = Chem.MolFromSmarts("[C;!R](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (fatty acid linkage) found"

    # Define a helper function that follows a strictly linear unbranched carbon chain.
    # It starts at a given carbon atom and returns the number of contiguous carbon atoms following the unique route.
    def linear_chain_length(atom, from_idx):
        # base length is 1 for the current atom
        length = 1
        # Find all neighboring carbon atoms; exclude the one we came from.
        nbr_carbons = [nbr for nbr in atom.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != from_idx]
        # If more than one neighbor carbon is found, this branch is not strictly linear.
        if len(nbr_carbons) != 1:
            return length
        next_atom = nbr_carbons[0]
        # For strict linearity, the next atom (if not terminal) must have exactly 2 carbon neighbors:
        # one is coming from current atom and one forward.
        # We then continue recursively.
        further_nbrs = [nbr for nbr in next_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
        if len(further_nbrs) > 1:
            return length + 1
        return length + linear_chain_length(next_atom, atom.GetIdx())

    # Set thresholds for an acceptable fatty acyl chain: minimum 12 and maximum 22 carbons.
    MIN_CARBONS = 12
    MAX_CARBONS = 22

    acyl_chain_found = False
    found_chain_length = 0

    # The thioester pattern gives a two-atom match: match[0] is the carbonyl carbon and match[1] the sulfur.
    for match in thioester_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Among neighbors of the carbonyl, choose a carbon that is not the sulfur.
        fatty_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_start = nbr
                break
        if fatty_start is None:
            continue  # no suitable carbon attached; try next match

        # Count the linear chain length starting from fatty_start.
        chain_length = linear_chain_length(fatty_start, carbonyl_idx)
        if MIN_CARBONS <= chain_length <= MAX_CARBONS:
            acyl_chain_found = True
            found_chain_length = chain_length
            break

    if not acyl_chain_found:
        return False, (
            f"Fatty acyl chain missing or not linear and within {MIN_CARBONS}-{MAX_CARBONS} carbons."
        )

    # Even if the overall formal charge isn't exactly -4 in the SMILES, we use our substructure evidence.
    return True, (
        f"Contains CoA moiety (adenine and pantetheine fragments), thioester group with a fatty acyl chain of "
        f"length {found_chain_length} carbons, and expected connectivity for long-chain fatty acyl-CoA(4-)."
    )

# Example usage:
if __name__ == "__main__":
    # Provide an example SMILES from the true positives
    test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)