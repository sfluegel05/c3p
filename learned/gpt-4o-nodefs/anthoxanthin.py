"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid characterized by two aromatic rings
    (often connected by a three-carbon bridge). They often have hydroxyl or methoxy 
    groups and can be glycosylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General flavonoid scaffolds could be like benzopyran structures
    flavonoid_patterns = [
        Chem.MolFromSmarts("c1cc2c(cc1)C=CO2"),  # Simple flavonoid scaffold
        Chem.MolFromSmarts("c1ccccc1-C2=C(O)c3ccccc3-OC2"),  # Generalized flavone structure
    ]
    
    # Check if any pattern matches
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_patterns):
        return False, "Does not contain typical anthoxanthin backbone structure"

    # Hydroxyl and methoxy groups are typical
    hydroxyl_or_methoxy_pattern = Chem.MolFromSmarts("[OX2H,OX1C]")
    if not mol.GetSubstructMatches(hydroxyl_or_methoxy_pattern):
        return False, "No hydroxyl or methoxy groups found, which are typical for anthoxanthins"

    # Optional: Check for possible glycoside linkages broadly
    glycoside_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](CO)O1")  # generic glycosidic pattern
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No evidence of glycosidic linkage that can appear in anthoxanthins, optional but common"

    return True, "Contains characteristic anthoxanthin backbone with potential hydroxyl/methoxy and glycosidic moieties"