"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group 
contains a double bond between positions 2 and 3.
"""

from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2‑enoyl‑CoA based on its SMILES string.
    A 2‑enoyl‑CoA is defined as an unsaturated fatty acyl-CoA in which the 
    acyl (fatty acid) moiety is attached via a thioester group and contains a double bond 
    between the alpha (position 2) and beta (position 3) carbons (i.e. between the carbonyl 
    and its adjacent carbon). In addition, the molecule should show features of CoA 
    (e.g. an adenine nucleotide).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 2‑enoyl‑CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, search for a thioester functional group using a SMARTS pattern.
    # The pattern "C(=O)S" will locate carbonyl groups directly attached to a sulfur.
    thioester_pat = Chem.MolFromSmarts("C(=O)S")
    if thioester_pat is None:
        return False, "Error building thioester SMARTS pattern"
    
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "Thioester group (C(=O)S) not found, so no acyl-CoA thioester"
    
    # For each thioester match, inspect the sulfur’s neighbor that is not the carbonyl carbon.
    # That neighbor is expected to be the acyl chain’s alpha carbon. We then check whether
    # that carbon participates in a carbon–carbon double bond (indicating a 2-enoyl group).
    unsat_acyl_found = False
    for match in thioester_matches:
        # 'match' is a tuple of indices corresponding to the SMARTS atom order:
        # index0: carbonyl carbon, index1: oxygen, index2: sulfur.
        s_idx = match[2]
        s_atom = mol.GetAtomWithIdx(s_idx)
        # Examine neighbors of the sulfur
        for nbr in s_atom.GetNeighbors():
            if nbr.GetIdx() == match[0]:
                # Skip the carbonyl carbon already part of the thioester
                continue
            # If the neighbor is a carbon, consider it as the acyl alpha carbon.
            if nbr.GetAtomicNum() == 6:
                # Loop over all bonds of this potential alpha carbon.
                for bond in nbr.GetBonds():
                    # Check if any bond is a double bond.
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # Make sure the double bond is to another carbon (the beta carbon)
                        other_atom = bond.GetOtherAtom(nbr)
                        if other_atom.GetAtomicNum() == 6:
                            unsat_acyl_found = True
                            break
                if unsat_acyl_found:
                    break
        if unsat_acyl_found:
            break

    if not unsat_acyl_found:
        return False, "The unsaturated acyl moiety (double bond between positions 2 and 3) was not found"
    
    # Next, confirm that the molecule contains features of CoA.
    # Here we look for the adenine ring typical for CoA structures.
    coa_pat = Chem.MolFromSmarts("n1cnc2ncnc12")
    if coa_pat is None:
        return False, "Error building CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pat):
        return False, "CoA structural features (adenine ring) not found in the molecule"
    
    return True, "Molecule contains a thioester-linked acyl group with a double bond between positions 2 and 3 and a CoA moiety"

# Example usage (uncomment to test):
# smi = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
# result, reason = is_2_enoyl_CoA(smi)
# print(result, reason)