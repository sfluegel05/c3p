"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is any molecule containing the isoflavonoid core structure,
    which is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a list of SMARTS patterns for isoflavonoid cores
    isoflavonoid_smarts_list = [
        # Isoflavone core (3-phenylchromen-4-one)
        'O=C1C=CC2=CC=CC=C2O1-c3ccccc3',
        # Isoflavanone core
        'O=C1CCC2=CC=CC=C2O1-c3ccccc3',
        # Isoflavene core
        'C1=CC=C2C(=C1)C=CC(O2)=C3C=CC=CC3',
        # Isoflavan core
        'O1CCC2=CC=CC=C2C1-c3ccccc3',
        # Pterocarpan core
        'O1CC[C@H]2Oc3c(OC2c1ccccc1)c(OC)c(OC)cc3C1=CC=CC=C1',
        # General benzopyran with aryl at position 3
        '[O]1C=CC=C2C=CC=CC12-c3cccc(c3)',
    ]

    # Search for isoflavonoid core structures
    for smarts in isoflavonoid_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # Skip invalid patterns
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains isoflavonoid core structure matching pattern: {smarts}"
    
    return False, "Does not contain isoflavonoid core structure"