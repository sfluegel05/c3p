from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxine(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxine.
    These are dibenzodioxines with 2 or more chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorinated dibenzodioxine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for dibenzodioxine core structure using SMARTS pattern
    # Two benzene rings connected by two oxygen bridges
    dioxin_pattern = Chem.MolFromSmarts('c1ccc2Oc3ccccc3Oc2c1')
    if not mol.HasSubstructMatch(dioxin_pattern):
        return False, "Does not contain dibenzodioxine core structure"
    
    # Count chlorine substituents
    chlorine_pattern = Chem.MolFromSmarts('[Cl]')
    chlorine_matches = mol.GetSubstructMatches(chlorine_pattern)
    chlorine_count = len(chlorine_matches)
    
    if chlorine_count < 2:
        return False, f"Only has {chlorine_count} chlorine substituents (minimum 2 required)"

    # Check that all non-oxygen, non-carbon atoms are chlorine
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'O', 'Cl']:
            return False, f"Contains non-allowed atom type: {atom.GetSymbol()}"

    # Check for exactly 2 oxygen atoms
    oxygen_count = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O'])
    if oxygen_count != 2:
        return False, f"Contains {oxygen_count} oxygen atoms (must be exactly 2)"

    return True, f"Polychlorinated dibenzodioxine with {chlorine_count} chlorine substituents"
# Pr=1.0
# Recall=1.0