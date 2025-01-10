"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    Biflavonoid typically contains two flavonoid moieties connected by a covalent bond.
    
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

    # Updated and broader pattern to search for flavonoid like connectivity
    # This tries to capture phenyl-benzopyran moieties connected in known ways in biflavonoids
    flavonoid_pattern = Chem.MolFromSmarts("c1ccccc1-c2c(O)cc3oc(=O)cc(=c3c2)-c2c(O)cc4oc(=O)cc(=c4c2)")

    # This ensures the presence of two flavonoid units connected typically by a single bond or through oxygen
    flavonoid_linker = Chem.MolFromSmarts("c1c2cc(O)cc3oc(=O)cc(=c3c2c(c1O)-c1ccccc1))")

    # Check if the molecule matches the flavonoid pattern twice and confirms suitable linkage
    match_basic = mol.HasSubstructMatch(flavonoid_pattern)
    has_linkage = mol.HasSubstructMatch(flavonoid_linker)

    if match_basic and has_linkage:
        return True, "Molecule contains dual connected flavonoid units typical of biflavonoids"
    
    return False, "Molecule does not match typical biflavonoid substructure patterns"