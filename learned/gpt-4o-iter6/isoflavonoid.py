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

    # Isoflavonoid pattern: 1-benzopyran structure with an aryl group at position 3
    isoflavonoid_pattern = Chem.MolFromSmarts("c1cc2ccc(c1)OCC2=c3cccc(c3)[#6]")

    if mol.HasSubstructMatch(isoflavonoid_pattern):
        return True, "Contains a 1-benzopyran core with an aryl group at position 3"
    
    # Check for additional features specific to isoflavonoids
    # Common functional groups like O-methylations, hydroxylations, etc.
    o_methoxy_pattern = Chem.MolFromSmarts("COc")
    if mol.HasSubstructMatch(o_methoxy_pattern):
        return True, "Contains isoflavonoid structure with methoxy group"

    hydroxyl_pattern = Chem.MolFromSmarts("O[cH]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return True, "Contains isoflavonoid structure with hydroxyl group"

    return False, "Does not contain the necessary isoflavonoid structure or functional groups"

# Example validation
examples = [
    "COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O",  # 2,3-dihydrobiochanin A
    "O=C1[C@H](COc2ccccc12)c1ccccc1"       # (3S)-isoflavanone
]

for e in examples:
    result, reason = is_isoflavonoid(e)
    print(f"SMILES: {e}, Result: {result}, Reason: {reason}")