"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python3
"""
Classifies: UDP-sugar
Definition:
  A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to
  an unspecified sugar via an anomeric diphosphate linkage.

This script uses RDKit to check for two main features:
  1) A pyrimidine nucleotide component: we require the presence of a uracil substructure.
  2) A diphosphate linkage: we check for an oxygen atom that bridges two phosphorus atoms 
     and is also connected to at least one carbon.
     
If both features are found, we classify the structure as a UDP-sugar.
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    A UDP-sugar must have:
      - A pyrimidine nucleotide component (here we require a uracil substructure), and
      - A diphosphate linkage where at least one oxygen atom bridges two phosphorus atoms 
        and also is connected to a carbon (suggesting attachment to a sugar).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is considered a UDP-sugar, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a uracil moiety (pyrimidine with two carbonyl groups)
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Uracil (nucleotide) substructure not found"
    
    # Look for a bridging oxygen in a diphosphate linkage.
    # The bridging oxygen should be attached to at least two phosphorus atoms and at least one carbon.
    diphosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:  # Only consider oxygen
            continue
        # Get neighboring phosphorus atoms
        p_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 15]
        if len(p_neighbors) < 2:
            continue  # Need at least two phosphorus neighbors to be bridging
        # Check if this oxygen is also attached to at least one carbon (attachment to a sugar)
        other_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if other_neighbors:
            diphosphate_found = True
            break

    if not diphosphate_found:
        return False, "Diphosphate linkage (bridging oxygen attached to two phosphorus and one carbon) not found"

    return True, "Molecule contains a uracil nucleotide and a diphosphate linkage attached to a sugar"

# Example usage: (you can remove or comment these out if integrating into a larger codebase)
if __name__ == "__main__":
    test_smiles = [
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)C(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",  # UDP-alpha-D-mannuronic acid
        "CC(=O)N[C@H]1[C@@H](O[C@H](CO)[C@@H](O)[C@@H]1OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(=O)[nH]c1=O"  # UDP-N-acetyl-beta-D-mannosaminouronic acid (example)
    ]
    
    for s in test_smiles:
        result, reason = is_UDP_sugar(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")