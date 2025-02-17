"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python3
"""
Classifies: UDP-sugar
Definition:
  A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to
  an unspecified sugar via an anomeric diphosphate linkage.
  
This script uses RDKit to check:
  1) That there is a pyrimidine nucleotide substructure (uracil).
  2) That there is a diphosphate linkage, which we now define by a P-O-P motif.
If both features are found, the molecule is classified as a UDP-sugar.
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    A UDP-sugar must have:
      - A pyrimidine nucleotide component (here defined as a uracil substructure), and
      - A diphosphate linkage (a P-O-P connectivity) that indicates a UDP moiety.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is considered a UDP-sugar, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a uracil moiety (pyrimidine with carbonyl groups).
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Uracil (nucleotide) substructure not found"
    
    # Define SMARTS for a diphosphate linkage.
    # This pattern matches a phosphorus atom connected via an oxygen to another phosphorus atom.
    diphosphate_smarts = "[P]-O-[P]"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate linkage (P-O-P) not found"

    return True, "Molecule contains a uracil substructure and a diphosphate (P-O-P) linkage indicative of a UDP-sugar"

# Example usage: (remove or comment these out if integrating into a larger codebase)
if __name__ == "__main__":
    test_smiles = [
        # UDP-alpha-D-mannuronic acid
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)C(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
        # A second UDP-sugar example (UDP-N-acetyl-beta-D-mannosaminouronic acid) might be added here.
    ]
    
    for s in test_smiles:
        result, reason = is_UDP_sugar(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")