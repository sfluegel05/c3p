"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: Diglyceride – a glyceride in which any two of the three hydroxyl groups are acylated.
Here we attempt to (1) identify a three‐carbon chain (a candidate glycerol backbone) that has exactly 
three oxygen substituents (one per carbon) and (2) classify the oxygen substituents as ester oxygen 
(via an acyl group attached, i.e. O–C(=O)–R) or as a free hydroxyl (an –OH group). 
If exactly two oxygens are esterified while one is free, we classify the structure as a diglyceride.
Note: This is a heuristic approach that may not catch every edge case.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    
    A diglyceride is defined as glycerol in which exactly two out of three hydroxy groups are acylated.
    This implementation tries to identify a glycerol backbone (a three-carbon linear chain with three 
    oxygen substituents) and then classifies the oxygen substituents as being either esterified (attached 
    to an acyl group via a carbonyl) or free hydroxyl.
    
    Args:
        smiles (str): SMILES string
        
    Returns:
        bool: True if the structure meets the diglyceride criteria, False otherwise.
        str: Reason or diagnostic message.
    """
    # Parse molecule; if invalid, return an error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject molecules with too low a molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a realistic diglyceride."
    
    # Add explicit hydrogens so that free hydroxyl groups will be seen.
    mol = Chem.AddHs(mol)
    
    # Gather candidate three-carbon chains (non-ring bonds only).
    atoms = list(mol.GetAtoms())
    candidate_found = False
    reason_details = ""
    
    # To reduce duplicates, we require a < b < c (by atom index).
    for a in atoms:
        if a.GetAtomicNum() != 6:
            continue
        a_idx = a.GetIdx()
        for b in a.GetNeighbors():
            if b.GetAtomicNum() != 6 or b.GetIdx() <= a_idx:
                continue
            # Check that bond a-b is not in a ring.
            bond_ab = mol.GetBondBetweenAtoms(a_idx, b.GetIdx())
            if bond_ab.IsInRing():
                continue
            for c in b.GetNeighbors():
                if c.GetAtomicNum() != 6 or c.GetIdx() <= b.GetIdx() or c.GetIdx() == a_idx:
                    continue
                bond_bc = mol.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx())
                if bond_bc.IsInRing():
                    continue

                # Candidate chain: atoms a, b, c (assuming linear order)
                chain_atoms = {a.GetIdx(), b.GetIdx(), c.GetIdx()}

                # For a glycerol-like backbone, we expect exactly one oxygen substituent on each carbon
                oxygen_substituents = {}  # key: glycerol carbon idx, value: the oxygen atom
                for atom in (a, b, c):
                    # Look at neighbors that are oxygen and NOT part of the chain.
                    oxy_neighbors = [nbr for nbr in atom.GetNeighbors() 
                                     if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in chain_atoms]
                    if len(oxy_neighbors) != 1:
                        # If not exactly one oxygen substituent, this chain is not a glycerol backbone.
                        break
                    oxygen_substituents[atom.GetIdx()] = oxy_neighbors[0]
                else:
                    # If we did not break out, then we have 3 oxygen substituents – one per carbon.
                    # Now classify each oxygen substituent as either ester oxygen or free hydroxyl.
                    ester_count = 0
                    free_oh_count = 0

                    for carbon in (a, b, c):
                        o_atom = oxygen_substituents[carbon.GetIdx()]
                        # For a free hydroxyl, we expect the oxygen to have an attached hydrogen.
                        # For an ester oxygen, the o_atom should be bound to a carbon that has a carbonyl.
                        o_neighbors = o_atom.GetNeighbors()
                        # Remove the glycerol carbon from its neighbor list.
                        attached = [nbr for nbr in o_neighbors if nbr.GetIdx() != carbon.GetIdx()]
                        
                        # Check for hydrogen neighbor.
                        has_H = any(nbr.GetAtomicNum() == 1 for nbr in o_neighbors)
                        if has_H:
                            free_oh_count += 1
                        else:
                            # If not free hydroxyl, then check if it is esterified:
                            # It should be attached to exactly one carbon (the acyl carbon).
                            if len(attached) == 1 and attached[0].GetAtomicNum() == 6:
                                acyl_carbon = attached[0]
                                # Check if acyl_carbon has a double bonded oxygen (i.e. a carbonyl). 
                                found_carbonyl = False
                                for nbr in acyl_carbon.GetNeighbors():
                                    bond = mol.GetBondBetweenAtoms(acyl_carbon.GetIdx(), nbr.GetIdx())
                                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        found_carbonyl = True
                                        break
                                if found_carbonyl:
                                    ester_count += 1
                                else:
                                    # If not, treat as free OH (this may catch some borderline cases)
                                    free_oh_count += 1
                            else:
                                # If the substitution at oxygen is ambiguous, count it as free.
                                free_oh_count += 1
                    # Now if the counts match exactly the diglyceride signature, we return True.
                    if ester_count == 2 and free_oh_count == 1:
                        return True, ("Matched candidate glycerol backbone (atoms %d-%d-%d) with "
                                      "esterified oxygens = 2 and free hydroxyl = 1." 
                                      % (a_idx, b.GetIdx(), c.GetIdx()))
                    else:
                        reason_details += ("Candidate backbone (atoms %d-%d-%d) had %d ester(s) and %d free OH(s). "
                                           % (a_idx, b.GetIdx(), c.GetIdx(), ester_count, free_oh_count))
    # If no candidate was accepted, include a diagnostic count using SMARTS.
    ester_pattern = Chem.MolFromSmarts("[O]-C(=O)")
    free_oh_pattern = Chem.MolFromSmarts("[OX2H]")
    num_esters = len(mol.GetSubstructMatches(ester_pattern)) if ester_pattern else 0
    num_free_oh = len(mol.GetSubstructMatches(free_oh_pattern)) if free_oh_pattern else 0
    diag = ("No valid glycerol backbone candidate found. " +
            f"Diagnostic: {num_esters} ester group(s) and {num_free_oh} free hydroxyl group(s) overall. " +
            reason_details)
    return False, diag

# Example usage for testing:
if __name__ == "__main__":
    test_examples = [
        # True positives – a few representative diglyceride examples.
        "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])O",
        "CCCCCCCCCCCCCCCC(=O)OC[C@H](CO)OC(=O)CCCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC",
        # One of the false negatives from before:
        "O(CCC([N+](C)(C)C)C([O-])=O)CC(O)COC(=O)CCCCCCC/C=C\\CCCCCCCC"
    ]
    for smi in test_examples:
        result, reason = is_diglyceride(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("----")