"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is defined as a nitrogen atom connected to three hydrocarbyl groups, which in SMILES terms means a nitrogen bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to obtain a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Simplified SMARTS pattern for tertiary amine: Nitrogen bonded to three carbon atoms
    # We ensure we are not excluding carbonyl or cyano group unnecessarily here, focusing on the N and its direct C connections.
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3](C)(C)C")
    
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Nitrogen bonded to three carbon atoms found indicating a tertiary amine"

    return False, "No nitrogen bonded to three carbon atoms found"

# Example asserts to test the function
assert is_tertiary_amine("CC(C)N(CC[C@H](c1ccccc1)c1cc(C)ccc1O)C(C)C")[0] == True  # tolterodine example
assert is_tertiary_amine("CCN(CC)CCOC(=O)C1(CCCCC1)C1CCCCC1")[0] == True  # dicyclomine example