"""
Classifies: CHEBI:37581 gamma-lactone
"""
#!/usr/bin/env python3
"""
Classifies: gamma-lactone
Definition: A lactone having a five-membered lactone ring.
A gamma-lactone is a cyclic ester formed by a five-membered ring made up 
of four carbons and one oxygen (the ring oxygen). One of the ring carbons 
(called the carbonyl carbon) bears a carbonyl group (C=O) where the O is 
external to the ring.
Improvement over prior attempt:
  • Instead of inspecting every five-membered ring manually, we use a SMARTS 
    pattern that captures the essential motif of a gamma-lactone: a 5‐membered 
    ring with the pattern C1(=O)CCC[O]1.
  • In this pattern, the ring consists of a carbonyl carbon (C1) substituted by a 
    double-bonded O (which is not considered “in‐ring”), three aliphatic carbon atoms,
    and one oxygen atom that is part of the ring.
Because no descriptor is perfect, there may still be edge cases, but this generally 
improves both the false positive and false negative rates.
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a molecule that contains a cyclic ester ring
    consisting of exactly five atoms (four carbons and one oxygen) with an exocyclic 
    carbonyl oxygen on the lactone carbon.
    
    This function uses a SMILES/SMARTS pattern ([#6]1(=O)[#6][#6][#6][O]1) that 
    essentially describes a five-membered ring with a carbonyl group attached to one 
    of the ring carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True if at least one gamma-lactone motif is found,
                     with an explanation message; otherwise False.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for gamma-lactone.
    # This pattern describes:
    #   a ring closure (1 ... 1) in which:
    #       - the first atom is a carbon (#6) that carries a double bond to
    #         an oxygen (=O) (this oxygen is NOT in the ring)
    #       - the next three atoms are any carbon (#6)
    #       - the last atom is the ring oxygen ([O])
    # Thus the ring is: C(=O)-C-C-C-O (with the =O exocyclic).
    gamma_lactone_smarts = "C1(=O)CCC[O]1"
    lactone_query = Chem.MolFromSmarts(gamma_lactone_smarts)
    if lactone_query is None:
        return False, "Error in SMARTS pattern"

    # Perform a substructure search to see if the molecule contains a gamma-lactone.
    matches = mol.GetSubstructMatches(lactone_query)
    if matches:
        # Found at least one match: now try to provide a meaningful explanation.
        return True, "Found a five‐membered lactone ring with a carbonyl group external to the ring."
    else:
        return False, "No gamma-lactone ring (five-membered cyclic ester with proper connectivity) found."
    
# Example usage (for testing purposes):
if __name__ == "__main__":
    # You can test one of the provided examples:
    test_smiles = "CCCCCC[C@H]1CCC(=O)O1"  # gamma-nonanolactone
    result, reason = is_gamma_lactone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)