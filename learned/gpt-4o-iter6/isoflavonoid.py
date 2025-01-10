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
        return (None, "Invalid SMILES string")

    # Define more comprehensive isoflavonoid patterns
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("c1cc2oc(C3=CC=CC=C3)cc2cc1"),        # General isoflavonoid core
        Chem.MolFromSmarts("c1cc2oc(C3=CC=C(O)C=C3)cc2cc1"),     # Hydroxy group on aryl substituent
        Chem.MolFromSmarts("c1cc2oc(C3=CC=C(OC)C=C3)cc2cc1"),    # Methoxy group on aryl substituent
        Chem.MolFromSmarts("c1cc2oc(C3=CC=C(OCC)C=C3)cc2cc1")    # Ether group on aryl substituent
    ]

    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches isoflavonoid substructure"

    # Consider more complex forms with sugar attachments
    glycosides_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O")
    if mol.HasSubstructMatch(glycosides_pattern):
        return True, "Contains isoflavonoid structure with glycoside attachments"

    return False, "Does not contain the necessary isoflavonoid structure or functional groups"

# Example validation
examples = [
    "COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O",  # 2,3-dihydrobiochanin A
    "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O"  # neobavaisoflavone
]

for e in examples:
    result, reason = is_isoflavonoid(e)
    print(f"SMILES: {e}, Result: {result}, Reason: {reason}")