"""
Classifies: CHEBI:17002 cholesteryl ester
"""
#!/usr/bin/env python
"""
Classifies: Cholesteryl esters
Definition: A sterol ester obtained by formal condensation of the carboxy group
of any carboxylic acid with the 3-hydroxy group of cholesterol.
This program uses a looser SMARTS pattern for the generic steroid (cholesterol)
nucleus with an oxygen substituent (which in cholesteryl esters is esterified)
and then verifies that the oxygen is bound to a carbonyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The method looks for a steroid (cholesterol-like) fused ring system with a substituent
    oxygen attached. In cholesteryl esters the free 3-hydroxy in cholesterol is esterified,
    meaning that oxygen is bound to an acyl group (i.e. a carbonyl, C=O).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is identified as a cholesteryl ester, False otherwise.
        str : Explanation of the classification decision.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a looser SMARTS pattern for a generic steroid nucleus.
    # The pattern tries to capture a fused tetracyclic ring system (three 6-membered and one 5-membered ring) 
    # with an oxygen substituent. We remove chiral specifications to be more permissive.
    chol_smarts = "C1CCC2CC(C3CCC4(C)CC(CC4)C3)CC2C1[O]"
    steroid_query = Chem.MolFromSmarts(chol_smarts)
    if steroid_query is None:
        return False, "Failed to parse steroid SMARTS pattern"
    
    # Look for substructure matches for the steroid nucleus with the oxygen substituent.
    matches = mol.GetSubstructMatches(steroid_query)
    if not matches:
        return False, "No steroid (cholesterol-like) nucleus with an oxygen substituent found"
    
    # For each match, we assume that the last atom in the match is the oxygen substituent.
    # Now check if that oxygen is involved in an ester linkage by being bound to a carbonyl group.
    for match in matches:
        o_idx = match[-1]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Examine each neighbor of this oxygen atom
        for nbr in o_atom.GetNeighbors():
            # We expect the oxygen to be attached to a carbon; note that the same oxygen might
            # have one bond to the steroid and one bond to the acyl carbon.
            if nbr.GetAtomicNum() == 6:
                # For the adjacent carbon, check for at least one double bond to an oxygen.
                for bond in nbr.GetBonds():
                    # Skip bond back to our oxygen atom from the steroid part.
                    other = bond.GetOtherAtom(nbr)
                    if other.GetIdx() == o_idx:
                        continue
                    # Check if the bond is a double bond and the other atom is oxygen.
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                        return True, ("Contains a steroid nucleus with an oxygen substituent that is part "
                                      "of an ester (esterified 3-hydroxy group) attached to an acyl carbonyl")
    
    return False, ("Steroid (cholesterol-like) nucleus found but no ester carbonyl attached to its "
                   "oxygen substituent; not a cholesteryl ester")

# Example usage
if __name__ == "__main__":
    # We test one example: cholesteryl linoleate.
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)[C@H](C)CCCC(C)C"
    result, reason = is_cholesteryl_ester(test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)