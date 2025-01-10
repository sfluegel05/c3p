"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid characterized by a heterocyclic pyrone 
    (or benzopyrone) backbone, often with several hydroxyl or methoxy groups and potential glycosylation.

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

    # Specific anthoxanthin flavonoid pattern: 1-benzopyran-4-one backbone
    anthoxanthin_pattern = Chem.MolFromSmarts("O=C1C=COC2=CC=CC=C12")  # Better representation of benzopyranone
    
    # Check if backbone structure matches
    if not mol.HasSubstructMatch(anthoxanthin_pattern):
        return False, "Does not contain typical anthoxanthin benzopyranone backbone structure"

    # Check for hydroxyl or methoxy groups
    hydroxyl_or_methoxy_pattern = Chem.MolFromSmarts("[OX2H,OX1C]")  # OH or OCH3 groups
    hydroxyl_or_methoxy_matches = mol.GetSubstructMatches(hydroxyl_or_methoxy_pattern)
    if len(hydroxyl_or_methoxy_matches) < 3:  # At least 3 such groups typically found in anthoxanthins
        return False, f"Insufficient hydroxyl or methoxy groups found, expected several typical of anthoxanthins"

    # Optional: Check for potential glycosidic linkages
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")  # Broad pattern
    if mol.HasSubstructMatch(glycoside_pattern):
        return True, "Anthoxanthin structure with potential glycosidic linkage"

    return True, "Contains characteristic anthoxanthin backbone and sufficient hydroxyl/methoxy groups"