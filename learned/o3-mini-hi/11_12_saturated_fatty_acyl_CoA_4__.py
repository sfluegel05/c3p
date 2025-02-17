"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
CHEBI:84947 is mentioned in the definition.
This classifier verifies that the molecule is a fatty acyl-CoA(4-) 
and that the fatty acyl portion has a saturated (i.e. single) bond between carbon 11 and 12.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-).
    
    A fatty acyl-CoA(4-) contains a fatty acyl chain attached via a thioester 
    (i.e. a [CX3](=O)[SX2] functional group) to a CoA moiety. 
    The 11-12-saturated fatty acyl-CoA(4-) class is defined as any such molecule in which
    the bond connecting the 11th and 12th carbon atoms (counting from the carbonyl carbon as C1)
    in the fatty acyl chain is saturated (a single bond).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): (True, reason) if molecule matches the class,
                     (False, reason) if not.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the CoA moiety.
    # We use a SMARTS pattern that is expected to exist in most acyl-CoA molecules.
    # This pattern looks for -SCCNC(=O)CCNC(=O), which is part of the CoA tail.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found - not a fatty acyl-CoA(4-) molecule"

    # Identify the thioester group: a carbonyl directly attached to a sulfur.
    # The SMARTS "[CX3](=O)[SX2]" should match the thioester linkage.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester group not found - molecule may not be fatty acyl-CoA(4-)"
    
    # We assume that the first thioester found is the fatty acyl thioester.
    # In the match, the first atom is the carbonyl carbon and the second is the sulfur.
    carbonyl_idx, S_idx = thioester_matches[0][0], thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # From the carbonyl, the fatty acyl chain is expected on the side that is not
    # the sulfur. Find the neighboring carbon (atomic num 6) that is not the S.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() != S_idx and nbr.GetAtomicNum() == 6:
            acyl_start = nbr
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found at the thioester group"

    # Build the linear fatty acyl chain.
    # According to fatty acid nomenclature, we count from the carbonyl carbon as C1.
    chain_atoms = [carbonyl_atom]  # C1
    previous_atom = carbonyl_atom
    current_atom = acyl_start         # C2
    
    # Traverse the chain until a branching occurs (or termination is reached).
    while True:
        chain_atoms.append(current_atom)
        # Look for a next carbon: choose neighbors that are carbons but not the previous.
        next_candidates = [nbr for nbr in current_atom.GetNeighbors()
                           if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous_atom.GetIdx()]
        # If exactly one next candidate, continue along a linear chain.
        if len(next_candidates) == 1:
            previous_atom, current_atom = current_atom, next_candidates[0]
        else:
            # End traversal at branching or termination.
            break

    # Check that the fatty acyl chain has at least 12 carbons (to contain an 11-12 bond).
    if len(chain_atoms) < 12:
        return False, f"Fatty acyl chain is too short (only {len(chain_atoms)} carbons) to have an 11-12 bond"

    # According to convention, numbering starts at the carbonyl carbon (C1).
    # The bond between the 11th and 12th carbon atoms will be between
    # chain_atoms[10] and chain_atoms[11] (0-indexed).
    atom_11 = chain_atoms[10]
    atom_12 = chain_atoms[11]
    bond_11_12 = mol.GetBondBetweenAtoms(atom_11.GetIdx(), atom_12.GetIdx())
    if bond_11_12 is None:
        return False, "11-12 bond not found in the fatty acyl chain"
    
    # Check the bond type: it should be a single bond (saturated) for this class.
    if bond_11_12.GetBondType() == Chem.BondType.SINGLE:
        return True, "The fatty acyl chain has a saturated (single) bond between carbon 11 and 12"
    else:
        return False, "The 11-12 bond is unsaturated (not a single bond), which does not meet the class criteria"
        
# Example usage:
if __name__ == "__main__":
    # Example SMILES from the provided list: (13Z)-3-oxodocosenoyl-CoA(4-)
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, ":", reason)