"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Flavonoids generally have a characteristic C6-C3-C6 structure, which may include various modifications.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broader flavonoid components with different oxidation states for the C ring
    patterns = [
        # General flavonoid pattern: A and B phenyl rings with C heterocyclic ring
        Chem.MolFromSmarts('[O=C]1c2ccccc2[o,c]3ccccc13'),  # Chromone core with flexible B-ring attachment
        Chem.MolFromSmarts('c1c2cc(O)c(O)cc2oc2ccccc12'),  # Flavonoid pattern with hydroxyl groups
        Chem.MolFromSmarts('[C@H]1(C=O)Cc2c(cc(O)c3OC1oc23)c1ccc(O)c(O)c1'),  # Flavanone/chromanone
        
        # Include common flavonoid glycosides recognized by high flexibility of the sugar unit
        Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O[C@@H]1O)-c2ccc(O)cc2'),  # Basic glycoside pattern

        # Recognize extensive flavonoid glycoside moieties
        Chem.MolFromSmarts('O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O-c2c3OCOc4cc(O)ccc4C(=O)c3ccc2'),  # Complex glycoside

        # Isoflavonoid backbone
        Chem.MolFromSmarts('Oc1ccc2c(c1)ccc(=O)o2'),  # Isoflavone

        # Include catechin patterns with open substitution accommodate
        Chem.MolFromSmarts('[O,c]1c(O)cc(O)[C@H](CC=O)[C@H]1c1ccc(O)c(O)c1'),  # Catechin/Anthocyanidin type
    ]

    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains flavonoid-like backbone or common modifications"

    return False, "No flavonoid-like backbone or pattern found"