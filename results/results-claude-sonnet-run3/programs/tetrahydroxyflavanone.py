from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydroxyflavanone(smiles: str):
    """
    Determines if a molecule is a tetrahydroxyflavanone (flavanone with at least 4 hydroxy groups).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tetrahydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for flavanone core structure
    # SMARTS pattern for flavanone core (C-C(=O)-C1-C=C-C=C-C1-O)
    flavanone_pattern = Chem.MolFromSmarts('O=C1CC(c2ccccc2)Oc2ccccc12')
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Not a flavanone - missing core structure"

    # Count number of hydroxy groups
    oh_pattern = Chem.MolFromSmarts('[OH]')
    hydroxy_matches = mol.GetSubstructMatches(oh_pattern)
    num_hydroxy = len(hydroxy_matches)

    # Check for minimum number of hydroxy groups
    if num_hydroxy < 4:
        return False, f"Only {num_hydroxy} hydroxy groups found - needs at least 4"

    return True, f"Tetrahydroxyflavanone with {num_hydroxy} hydroxy groups"
# Pr=0.8571428571428571
# Recall=1.0