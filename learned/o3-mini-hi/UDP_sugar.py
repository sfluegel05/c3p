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
  1) A uracil moiety (as found in the UDP nucleotide part).
  2) A diphosphate linkage where an oxygen (that bridges two phosphorus atoms)
     is also linked to a carbon (suggesting attachment to a sugar).
     
If both features are found, we classify the structure as a UDP-sugar.
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    A UDP-sugar must have:
      - a pyrimidine nucleotide component (here we require a uracil substructure), and
      - a diphosphate linkage where an oxygen bridges two phosphorus atoms and
        is connected to a carbon (suggesting attachment to a sugar).
    
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
    
    # SMARTS pattern for a uracil moiety
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Uracil (nucleotide) substructure not found"
    
    # Now search for a diphosphate linkage.
    #
    # We will search for an oxygen atom that connects two phosphorus atoms.
    # Then we check that the bridging oxygen is also attached (besides the P atoms)
    # to a carbon (suggesting a bond to an organic sugar group).
    diphosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen neighbor of phosphorus
                    # Find out if this oxygen is also connected to another phosphorus.
                    connected_ps = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 15 and n.GetIdx() != atom.GetIdx()]
                    if connected_ps:
                        # Now check that this bridging oxygen (nbr) is also attached (aside from the 2 phosphorus neighbors)
                        # to at least one carbon atom (atomic num 6).
                        other_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() not in (15,)]
                        if any(n.GetAtomicNum() == 6 for n in other_neighbors):
                            diphosphate_found = True
                            break
            if diphosphate_found:
                break
    if not diphosphate_found:
        return False, "Diphosphate linkage (with bridging O attached to a carbon) not found"
    
    # Optional extra check: sometimes one might want to verify that a cyclic sugar (a ring with oxygen atoms)
    # is present. This is more challenging because sugars vary widely. For now we assume that if both
    # a uracil and the proper diphosphate linkage are present, the structure is most likely a UDP-sugar.
    
    return True, "Molecule contains a uracil nucleotide and a diphosphate linkage attached to a sugar"

# Example usage: (you can remove or comment these out if integrating into a larger codebase)
if __name__ == "__main__":
    # A few example SMILES from the UDP-sugar class (taken from the user specification)
    test_smiles = [
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)C(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",  # UDP-alpha-D-mannuronic acid
        "CC(=O)N[C@H]1[C@@H](O[C@H](CO)[C@@H](O)[C@@H]1OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(=O)[nH]c1=O"  # UDP-N-acetyl-beta-D-mannosaminouronic acid (example)
    ]
    
    for s in test_smiles:
        result, reason = is_UDP_sugar(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")