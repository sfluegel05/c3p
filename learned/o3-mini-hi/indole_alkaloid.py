"""
Classifies: CHEBI:38958 indole alkaloid
"""
#!/usr/bin/env python
"""
Classifies: Indole Alkaloid - An alkaloid containing an indole skeleton.
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is defined as an alkaloid containing an indole skeleton,
    which is the bicyclic combination of a benzene ring fused to a pyrrole ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains an indole skeleton, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the indole scaffold.
    # This pattern captures a pyrrole ring (with the ring nitrogen) fused to a benzene ring.
    indole_pattern = Chem.MolFromSmarts("n1c2ccccc2c1")
    
    # Check if the molecule contains the indole substructure.
    if mol.HasSubstructMatch(indole_pattern):
        return True, "Molecule contains an indole skeleton characteristic of indole alkaloids."
    else:
        return False, "Molecule does not contain an indole skeleton."
        
# Example usage (uncomment for testing):
# test_smiles = "Oc1ccc2c(c1)c1ccnc3ccc(=O)n2c13"  # 10-hydroxycanthin-6-one
# result, reason = is_indole_alkaloid(test_smiles)
# print(result, reason)