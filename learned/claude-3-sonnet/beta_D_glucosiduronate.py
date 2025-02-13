"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:27718 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a carbohydrate acid derivative anion obtained by deprotonation
    of the carboxy group of any beta-D-glucosiduronic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for deprotonated carboxyl group(s)
    deprotonated_pattern = Chem.MolFromSmarts("[O-,C-](=O)O")
    deprotonated_matches = mol.GetSubstructMatches(deprotonated_pattern)
    if not deprotonated_matches:
        return False, "No deprotonated carboxyl group found"

    # Look for glucuronic acid moiety
    glucuronic_acid_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)C(=O)[O-])O)O)O"),  # Standard pattern
        Chem.MolFromSmarts("O[C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)C(=O)O)O)O)O"),  # Alternate pattern
        # Add more patterns as needed
    ]
    glucuronic_acid_match = False
    for pattern in glucuronic_acid_patterns:
        if mol.HasSubstructMatch(pattern):
            glucuronic_acid_match = True
            break
    if not glucuronic_acid_match:
        return False, "No glucuronic acid moiety found"

    # Check for additional structural features (optional)
    # ...

    return True, "Contains a deprotonated carboxyl group and a glucuronic acid moiety"