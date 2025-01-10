"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Biflavonoids typically consist of two flavonoid units. We'll define a simple pattern to look for such units.
    # This is a general check and may not match all biflavonoids precisely.
    flavonoid_pattern_1 = Chem.MolFromSmarts("c1cc(O)c2c(c1)oc(=O)cc2-c1c(O)cc(O)cc1")  # Example flavonoid unit
    flavonoid_pattern_2 = Chem.MolFromSmarts("Oc1cccc2c1c(=O)c(O)cc2-c1c(O)cc(O)cc1")  # Example flavonoid unit

    # Check if the molecule matches the flavonoid pattern at least twice
    match_1 = mol.HasSubstructMatch(flavonoid_pattern_1)
    match_2 = mol.HasSubstructMatch(flavonoid_pattern_2)

    if match_1 and match_2:
        return True, "Molecule contains two flavonoid units typical of biflavonoids"
    
    return False, "Molecule does not match typical biflavonoid substructure patterns"