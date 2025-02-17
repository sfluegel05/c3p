"""
Classifies: CHEBI:37739 glycerophospholipid
"""
#!/usr/bin/env python
"""
Classifies: glycerophospholipid
Definition:
  Any glycerolipid having a phosphate group ester-linked to a terminal (primary) carbon 
  of the glycerol backbone (i.e. one that is CH2 and only attached to one other carbon)
  and at least one fatty acyl chain connected as an acyl ester (OC(=O)) group.
Heuristic criteria used in this revision:
  1. Locate a phosphate group (P) that is connected via an oxygen (O) to a carbon which is “terminal”:
     • that carbon must have exactly two hydrogens (CH2) and only one neighbor that is a carbon,
       to increase the likelihood that it is the terminal carbon of a glycerol motif.
  2. Ensure at least one acyl ester group (OC(=O)) is present that is not part of the phosphate group.
  3. Relax the molecular weight filter (here ≥400 Da) to better capture lower‐mass glycerophospholipids.
Note: Because of the chemical diversity in glycerophospholipid structures, this is a heuristic filter.
It may still give some false positives and negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is defined (heuristically) as a glycerolipid having a phosphate 
    group ester-linked to a terminal (primary, CH2) carbon of the glycerol backbone plus at least
    one acyl ester (OC(=O)) group indicating a fatty-acyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a glycerophospholipid; False otherwise.
        str: Explanation/reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Look for a phosphate group (P) with an attached oxygen that is connected to a terminal carbon.
    phosphate_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check for a P-O bond (in either direction)
        if (a1.GetSymbol() == 'P' and a2.GetSymbol() == 'O') or (a1.GetSymbol() == 'O' and a2.GetSymbol() == 'P'):
            # Identify the oxygen atom from the bond
            o_atom = a1 if a1.GetSymbol() == 'O' else a2
            # Ensure that this oxygen is not the one that is double-bonded to P (if any)
            # Now, look for a carbon neighbor (other than P) attached to this oxygen.
            for nbr in o_atom.GetNeighbors():
                if nbr.GetSymbol() == 'C':
                    # Check if the carbon has exactly 2 implicit/explicit hydrogens (a CH2),
                    # and that it is terminal (only one C neighbor).
                    if nbr.GetTotalNumHs() == 2:
                        c_neighbors = [nb for nb in nbr.GetNeighbors() if nb.GetSymbol() == 'C']
                        if len(c_neighbors) == 1:
                            phosphate_found = True
                            break
            if phosphate_found:
                break
    if not phosphate_found:
        return False, "No phosphate group found that is ester-linked to a terminal CH2 carbon as expected in a glycerol backbone"
    
    # Criterion 2: Look for at least one acyl ester group (OC(=O)) that is not part of the phosphate.
    # We define the ester pattern by a SMARTS: an oxygen singly bonded to a carbonyl carbon.
    # (The [!P] filter helps avoid matching phosphate-connected oxygens.)
    ester_smarts = "[O;!$(*~P)]C(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No acyl ester (OC(=O)) group found that would indicate a fatty-acyl chain"
    
    # Criterion 3: Check molecular weight, but relax compared to previous threshold;
    # many glycerophospholipids are lower than 500 Da. Here we use 400 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight is too low to be a glycerophospholipid"
    
    # If all the criteria are met, classify as glycerophospholipid.
    return True, "Molecule contains a phosphate group ester-linked to a terminal CH2 of a glycerol backbone and an acyl ester group, consistent with a glycerophospholipid structure"

# Example usage (for testing):
if __name__ == "__main__":
    # Test one of the provided glycerophospholipid SMILES strings.
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC(=O)\\C=C\\C(O)=O"
    result, reason = is_glycerophospholipid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)