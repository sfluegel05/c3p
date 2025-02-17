"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty acyl-CoA(4-)
Definition: A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups
of any long-chain fatty acyl-CoA; major species at pH 7.3.
This implementation verifies that the molecule contains (a) an adenine nucleotide fragment (using several SMARTS variants),
(b) a pantetheine fragment, (c) a thioester group linking a fatty acyl chain and (d) a linear (unbranched)
hydrocarbon chain of 12–22 carbons attached to the thioester.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines whether a given SMILES string belongs to the class long-chain fatty acyl-CoA(4-).

    Process:
      1. Parse the SMILES string.
      2. Verify that an adenine nucleotide fragment is present (using several alternative SMARTS).
      3. Check for a pantetheine fragment characteristic of CoA.
      4. Find a thioester group ([C;!R](=O)[S]) indicating the fatty acyl linkage.
      5. From the thioester carbonyl carbon, trace the attached carbon (not the one bound to S) into a chain.
         • The chain-ending routine is iterative and will only continue if a unique carbon neighbor is found
           (enforcing linear—unbranched—connectivity).
      6. Accept if at least one thioester leads into a linear carbon chain of 12–22 carbons.

    Args:
        smiles (str): SMILES string of the candidate molecule.

    Returns:
        bool: True if the molecule meets all criteria; False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Step 1: Check for adenine nucleotide fragment ---
    # Use several SMARTS patterns to allow slightly different representations.
    adenine_smarts = [
        "n1cnc2ncnc12",         # plain adenine core
        "n1cnc2c(n)ncnc12",      # variant with a methyl substituent
        "n1cnc2c(N)ncnc12"       # variant with an NH2 substituent on one ring
    ]
    found_adenine = False
    for smarts in adenine_smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            found_adenine = True
            break
    if not found_adenine:
        return False, "Adenine nucleotide fragment not found; not a CoA derivative"

    # --- Step 2: Check for the pantetheine fragment ---
    # This SMARTS should match the key connectivity of pantetheine.
    pantetheine_patt = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not pantetheine_patt or not mol.HasSubstructMatch(pantetheine_patt):
        return False, "Pantetheine fragment not found; not a CoA derivative"

    # --- Step 3: Locate a thioester group ---
    # A thioester is defined as a non-ring carbon (C;!R) with a =O and a single bond to a sulfur.
    thioester_patt = Chem.MolFromSmarts("[C;!R](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_patt)
    if not thioester_matches:
        return False, "No thioester group (fatty acid linkage) found"

    # --- Step 4: Once a thioester is found, follow the fatty acyl chain ---
    # This iterative function walks from a given carbon atom into a chain.
    def trace_linear_chain(start_atom, from_atom):
        """Trace a chain iteratively starting from start_atom;
           at each step only continue if there is exactly one carbon neighbor (excluding the atom we came from).
           Returns the chain length (number of carbon atoms encountered, including start_atom)."""
        length = 1  # start_atom counts as one carbon
        prev_atom = from_atom
        current_atom = start_atom
        while True:
            # Get neighboring carbons except the one we came from.
            nbr_carbons = [nbr for nbr in current_atom.GetNeighbors() 
                           if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()]
            # For a linear chain exactly one such neighbor must exist.
            if len(nbr_carbons) != 1:
                break
            # If the current atom itself has extra carbon joiners (i.e. a branch), do not continue.
            # (We count only if the continuation is unique.)
            next_atom = nbr_carbons[0]
            length += 1
            prev_atom = current_atom
            current_atom = next_atom
        return length

    MIN_CHAIN = 12
    MAX_CHAIN = 22
    acyl_chain_valid = False
    chain_length_found = 0

    # For each thioester, try to find a fatty acyl chain.
    for match in thioester_matches:
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From the carbonyl, choose a neighbor that is a carbon and not the sulfur.
        fatty_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_start = nbr
                break
        if fatty_start is None:
            continue  # no candidate fatty acyl carbon found with this thioester; try next

        # Trace the potential fatty chain.
        chain_length = trace_linear_chain(fatty_start, carbonyl_atom)
        if MIN_CHAIN <= chain_length <= MAX_CHAIN:
            acyl_chain_valid = True
            chain_length_found = chain_length
            break

    if not acyl_chain_valid:
        return False, f"Fatty acyl chain missing or not linear with required {MIN_CHAIN}-{MAX_CHAIN} carbons."

    # If we reach here then the molecule meets all criteria.
    return True, (
        f"Contains CoA moiety (adenine and pantetheine fragments), a thioester group with a fatty acyl chain of length {chain_length_found} carbons, "
        "and expected connectivity for long-chain fatty acyl-CoA(4-)."
    )


# For testing purposes:
if __name__ == "__main__":
    # Example SMILES from one of the true positive cases.
    test_smiles = "CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, explanation = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, explanation)