"""
Classifies: CHEBI:17002 cholesteryl ester
"""
#!/usr/bin/env python
"""
Classifies: Cholesteryl esters
Definition: A sterol ester obtained by formal condensation of the carboxy group
of any carboxylic acid with the 3-hydroxy group of cholesterol.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The method first looks for a cholesterol substructure
    (i.e. the steroid fused ring system with an oxygen substituent at the 3-position).
    Then it checks that the oxygen is involved in an ester linkage by being bound to
    a carbon that is double-bonded to an oxygen (i.e. a carbonyl group).
    
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
    
    # Define a SMARTS pattern for the cholesterol skeleton.
    # The pattern is based on the cholesterol core with the 3-hydroxy group.
    # In cholesteryl esters the hydroxyl oxygen becomes esterified (losing its H)
    # but the oxygen remains present. This SMARTS pattern matches that oxygen.
    chol_smarts = ("C1CC[C@H]2CC[C@H]3CC[C@]4(C)[C@@H](CC[C@H]4[C@@H]3CC[C@]12C)[O]")
    chol_query = Chem.MolFromSmarts(chol_smarts)
    if chol_query is None:
        return False, "Failed to parse cholesterol SMARTS pattern"
    
    # Look for substructure matches for the cholesterol moiety.
    matches = mol.GetSubstructMatches(chol_query)
    if not matches:
        return False, "No cholesterol skeleton with the expected 3-hydroxy (esterified) group found"

    # For each match, the last atom in the SMARTS match (the [O]) should be the ester oxygen.
    # Check if this oxygen is attached to a carbonyl group.
    for match in matches:
        # The oxygen atom from the cholesterol part is captured as the last atom in the SMARTS pattern.
        o_idx = match[-1]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Go over each neighbor of this oxygen.
        for nbr in o_atom.GetNeighbors():
            # We expect the oxygen (of the ester linkage) to be bonded to a carbon
            # that is part of the acyl (fatty acid) moiety.
            if nbr.GetAtomicNum() == 6:  # Check for carbon
                # For this neighboring carbon, check if it is part of a carbonyl group.
                for bond in nbr.GetBonds():
                    # Skip bond back to the original oxygen
                    other = bond.GetOtherAtom(nbr)
                    if other.GetIdx() == o_idx:
                        continue
                    # Check if the bond is a double bond and the other atom is oxygen.
                    # This is our indicator of a carbonyl.
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                        return True, "Contains a cholesterol moiety with its 3-hydroxy group esterified to an acyl (carbonyl) group"
    
    return False, "Cholesterol skeleton found but no ester carbonyl attached to the 3-hydroxy oxygen"

# Example usage:
if __name__ == "__main__":
    # Example: cholesteryl linoleate
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)[C@H](C)CCCC(C)C"
    result, reason = is_cholesteryl_ester(test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)