"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: long-chain fatty acyl-CoA(4-)
Definition: A fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups 
of any long-chain fatty acyl-CoA; major species at pH 7.3.
This program uses simple substructure searches to detect a CoA moiety and a thioester group from which
the acyl chain is identified. The fatty acyl chain is assumed to be “long‐chain” if it contains at least 
12 contiguous carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) belongs to the class
    long-chain fatty acyl-CoA(4-).

    We check for three features:
      1. The presence of a CoA moiety (using a characteristic fragment of CoA).
      2. A thioester group (C(=O)S) which links the fatty acyl chain with the CoA.
      3. A long fatty acyl chain: we count the number of contiguous carbons attached 
         to the carbonyl carbon (threshold: at least 12 carbons).
      4. Overall formal charge of -4 (as expected for CoA(4-)).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as long-chain fatty acyl-CoA(4-), False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check overall formal charge (expecting -4 for the deprotonated CoA species)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Formal charge is {total_charge} (expected -4 for CoA(4-))"
    
    # Check for the CoA moiety using a recognizable fragment.
    # Here we use a fragment pattern from the pantetheine and nucleotide parts: "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Look for a thioester group: pattern for a carbonyl carbon (with degree 3) bound to a sulfur.
    # The SMARTS "[C;D3](=O)[S]" should capture the thioester group.
    thioester_pattern = Chem.MolFromSmarts("[C;D3](=O)[S]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (fatty acid linkage) found"

    # Define a recursive function to count the length of a contiguous carbon chain.
    # We only travel along carbon atoms and do not branch (if branching is encountered, we take the longest route).
    def count_chain_length(atom, from_idx, visited):
        # Count current atom (note: visited is used to avoid loops)
        length = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == from_idx or nbr.GetIdx() in visited:
                continue
            # We only continue if the neighbor is a carbon
            if nbr.GetAtomicNum() == 6:
                visited.add(nbr.GetIdx())
                branch_length = 1 + count_chain_length(nbr, atom.GetIdx(), visited)
                if branch_length > length:
                    length = branch_length
        return length

    # Try each thioester match to find a fatty acyl chain
    # In the pattern "[C;D3](=O)[S]", match[0] is the carbonyl carbon and match[1] is the sulfur.
    chain_found = False
    acyl_chain_length = 0
    for match in thioester_matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Look among the neighbors of the carbonyl carbon for a carbon that is not the S from the thioester.
        fatty_start = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != match[1]:
                fatty_start = nbr
                break
        if fatty_start is None:
            continue  # no attached carbon found; try next match
        
        # Count the acyl chain length. We start from fatty_start.
        # We assume the chain is linear so we take the maximal contiguous carbon count.
        chain_length = count_chain_length(fatty_start, carbonyl_idx, {fatty_start.GetIdx()})
        if chain_length >= 12:
            chain_found = True
            acyl_chain_length = chain_length
            break

    if not chain_found:
        return False, "Fatty acyl chain not long enough (requires at least 12 contiguous carbons)"

    return True, (
        f"Contains CoA moiety, thioester group with a fatty acyl chain of length {acyl_chain_length} carbons, "
        "and overall charge -4"
    )

# Example usage:
if __name__ == "__main__":
    test_smiles = "CCCC\\C=C/C\\C=C/CCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_long_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)