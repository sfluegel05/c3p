"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Isoflavonoids
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at the 3-position.

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

    # Look for the isoflavonoid pattern (benzopyran with aryl at position 3)
    # Isoflavonoids are 1-benzopyran derivatives with an aryl group at the 3 position.
    
    # Isoflavonoid SMARTS pattern
    isoflavonoid_pattern = Chem.MolFromSmarts("c1cc2oc(cc3ccc[nH]3)c3ccccc3c2c1")  # Adjust to better capture isoflavonoid structure
    
    if mol.HasSubstructMatch(isoflavonoid_pattern):
        return True, "Contains a 1-benzopyran core with an aryl group at position 3"
    else:
        return False, "Does not contain the necessary isoflavonoid structure"

# Example validation
examples = [
    "COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O",  # 2,3-dihydrobiochanin A
    "O=C1[C@H](COc2ccccc12)c1ccccc1"       # (3S)-isoflavanone
]

for e in examples:
    result, reason = is_isoflavonoid(e)
    print(f"SMILES: {e}, Result: {result}, Reason: {reason}")