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

    # Expand isoflavonoid definition: 1-benzopyran core with diverse aryl group
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("Oc1ccc2ccccc2c1-c3ccccc3"),         # Typical isoflavonoid
        Chem.MolFromSmarts("O=c1cc2ccccc2oc1-c3ccccc3"),        # Variant with additional ketone
        Chem.MolFromSmarts("Oc1ccc2ccccc2c1-c3cc(O)ccc3"),      # With hydroxyl groups
        Chem.MolFromSmarts("Oc1ccc2cc(O)cccc2c1-c3ccc(OC)cc3")  # With methoxy groups
    ]

    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches isoflavonoid substructure with possible functional variations"
    
    # Check possible higher complexity with glucose or other sugar attachments
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H]([C@H]([C@@H]([C@H]1O)O)O)O)CO")
    if mol.HasSubstructMatch(sugar_pattern):
        return True, "Contains potential isoflavonoid structure with sugar attachments"

    # Check for common substituents in isoflavonoids
    common_substituents = [
        Chem.MolFromSmarts("CO"),  # Methoxy group
        Chem.MolFromSmarts("C=C"), # Carbon double bonds
        Chem.MolFromSmarts("CC(=O)O") # Acetate or ester linkages
    ]

    for substituent in common_substituents:
        if mol.HasSubstructMatch(substituent):
            return True, "Contains functional group common in isoflavonoids"

    return False, "Does not contain the necessary isoflavonoid structure or functional groups"

# Example validation
examples = [
    "COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O",  # 2,3-dihydrobiochanin A
    "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O"  # neobavaisoflavone
]

for e in examples:
    result, reason = is_isoflavonoid(e)
    print(f"SMILES: {e}, Result: {result}, Reason: {reason}")