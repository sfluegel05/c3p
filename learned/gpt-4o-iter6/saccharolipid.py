"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is expected to contain both a carbohydrate moiety and a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string provided."

    # Focused patterns for carbohydrate moieties typical in saccharolipids
    carbohydrate_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),  # Classic sugar rings
        Chem.MolFromSmarts("O[C@H]1O[C@H](C(O)C1O)C(O)=O"),              # Modified sugars (e.g., acidic sugars)
        Chem.MolFromSmarts("O[C@H]([C@H](O)[C@H](CO)O)CO")               # Sugar with free hydroxyl groups
    ]

    # Focused pattern for lipid components that can be considered distinct from generic glycolipids
    lipid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O"),    # Long chain fatty acid ester
        Chem.MolFromSmarts("CCCCCCCCCCCC(=O)N"),      # Long chain fatty acid amide
        Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCC"),     # Long saturated hydrocarbon chain
    ]

    # Check for a carbohydrate component
    has_carbohydrate = any(mol.HasSubstructMatch(carb_pattern) for carb_pattern in carbohydrate_patterns)
    if not has_carbohydrate:
        return False, "No characteristic carbohydrate moiety found."

    # Check for a lipid component
    has_lipid = any(mol.HasSubstructMatch(lipid_pattern) for lipid_pattern in lipid_patterns)
    if not has_lipid:
        return False, "No suitable lipid component found."

    # Ensure the carbohydrate and lipid are distinct parts
    for lipid_pattern in lipid_patterns:
        for carb_pattern in carbohydrate_patterns:
            if mol.HasSubstructMatch(lipid_pattern) and mol.HasSubstructMatch(carb_pattern):
                # Ensure these are distinct structural parts
                lipid_matches = mol.GetSubstructMatches(lipid_pattern)
                carb_matches = mol.GetSubstructMatches(carb_pattern)
                overlap = any(set(lm).intersection(set(cm)) for lm in lipid_matches for cm in carb_matches)
                if overlap:
                    continue

                return True, "Contains both distinct carbohydrate and lipid components, classifiable as a saccharolipid."

    return False, "Overlapping detected structures; does not fulfill distinct saccharolipid criteria."