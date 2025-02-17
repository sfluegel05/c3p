"""
Classifies: CHEBI:22315 alkaloid
"""
#!/usr/bin/env python3
"""
Classifies: Alkaloid

Definition:
    Alkaloids are naturally occurring, basic nitrogen compounds,
    usually containing one or more heterocyclic rings. They are
    primarily found in the plant kingdom (and also in bacteria,
    fungi, and animals). In contrast, compounds in which the nitrogen
    is exclusively exocyclic (such as dopamine, mescaline, serotonin, etc.)
    are generally classed as amines rather than alkaloids.
    
This function uses a simple heuristic:
  1. It verifies that the input SMILES can be converted to an RDKit molecule.
  2. It checks that the molecule contains at least one nitrogen atom.
  3. It inspects all nitrogen atoms to see if at least one is part of a ring.
If a nitrogen is found as part of a ring then the molecule is flagged as a potential alkaloid.
Otherwise, it is not classified as an alkaloid.
"""

from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string using a heuristic.
    
    An alkaloid (by our heuristic) is defined as a naturally occurring, basic nitrogen compound 
    that contains at least one nitrogen atom incorporated into a ring system (heterocycle).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alkaloid, False otherwise.
        str: Explanation of the classification reasoning.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all atoms that are nitrogen (atomic number 7)
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found; alkaloids require at least one nitrogen"

    # Obtain ring information from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Check if at least one nitrogen atom is part of a ring
    nitrogen_in_ring = False
    for atom in nitrogen_atoms:
        idx = atom.GetIdx()
        # Check each ring to see if the current nitrogen atom is part of it
        for ring in ring_info:
            if idx in ring:
                nitrogen_in_ring = True
                break
        if nitrogen_in_ring:
            break

    if not nitrogen_in_ring:
        return False, "All nitrogen atoms are exocyclic; likely an amine, not an alkaloid"
    
    # We could add more advanced checks such as ensuring basicity, molecular weight, or excluding peptides,
    # but in this heuristic presence of an in‐ring nitrogen is a strong indicator for many alkaloids.
    
    return True, "Contains a nitrogen incorporated in a ring system; matches key features of alkaloids"
    
# For testing purposes, you might run:
if __name__ == "__main__":
    test_smiles = "C[C@H](Cc1ccccc1)N(C)CC#C"  # Example (-)-selegiline; note this one has an exocyclic N.
    result, reason = is_alkaloid(test_smiles)
    print("Result for test molecule:", result)
    print("Reason:", reason)