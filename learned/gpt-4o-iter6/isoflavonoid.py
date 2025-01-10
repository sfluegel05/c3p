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
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

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

    # Define SMARTS patterns for isoflavonoid cores and substitutions
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("c1cc2oc(c3ccccc3)cc2cc1"),          # Basic isoflavonoid core structure
        Chem.MolFromSmarts("c1cc2oc(c3ccc(cc3)O)cc2cc1"),       # Hydroxy substitution on aryl group
        Chem.MolFromSmarts("c1cc2oc(c3ccc(cc3)OC)cc2cc1"),      # Methoxy substitution on aryl group
        Chem.MolFromSmarts("[O;H1]c1ccc2c(c1)cc(c3ccccc3)o2"),  # Free phenolic hydroxyl in 3-position
        Chem.MolFromSmarts("C1COc2cc(ccc2O)C=1"),               # Isoflavan structure
    ]

    glycosides_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O")
    
    # Check for isoflavonoid basic and substituted structures
    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            if mol.HasSubstructMatch(glycosides_pattern):
                return True, "Contains isoflavonoid core with glycoside attachments"
            return True, "Matches isoflavonoid core structure"

    return False, "Does not contain the necessary isoflavonoid structure or functional groups"

# Example validation
examples = [
    "COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O",  # 2,3-dihydrobiochanin A
    "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O"  # neobavaisoflavone
]

for e in examples:
    result, reason = is_isoflavonoid(e)
    print(f"SMILES: {e}, Result: {result}, Reason: {reason}")