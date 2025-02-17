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
    An indole alkaloid is defined as an alkaloid containing an indole skeleton:
    a bicyclic system consisting of a benzene ring fused with a pyrrole ring.
    
    Our approach first parses the SMILES string to a molecule and then uses a SMARTS pattern
    that captures the indole scaffold. We use the SMARTS pattern "c1ccc2c(c1)[nX2]cc2" which
    matches an aromatic six‐membered ring fused to a five‐membered ring that contains an aromatic
    nitrogen which may be substituted.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an indole skeleton, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an indole scaffold.
    # "c1ccc2c(c1)[nX2]cc2" matches a benzene ring fused to a five-membered heterocycle
    # containing an aromatic nitrogen (which may be substituted).
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nX2]cc2")
    if indole_pattern is None:
        return False, "Failed to create SMARTS pattern"
    
    # Check if the molecule contains the indole substructure.
    if mol.HasSubstructMatch(indole_pattern):
        return True, "Molecule contains an indole skeleton characteristic of indole alkaloids."
    else:
        return False, "Molecule does not contain an indole skeleton."

# Example usage (for testing, uncomment the following lines):
# test_smiles = "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)N1C3=C(C=CC=C3)C3=C1C1=C(C4=C(C=CC=C4)N21)C1=C3CNC1=O"  # staurosporine
# result, reason = is_indole_alkaloid(test_smiles)
# print(result, reason)