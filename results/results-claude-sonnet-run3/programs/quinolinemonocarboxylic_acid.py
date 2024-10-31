from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_quinolinemonocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a quinoline monocarboxylic acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinoline monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    if len(carboxyl_matches) > 1:
        return False, "More than one carboxylic acid group found"

    # Check for quinoline core structure
    # Quinoline is a bicyclic structure with a benzene ring fused to a pyridine ring
    quinoline_pattern = Chem.MolFromSmarts('c1cccc2c1[n]ccc2')  # Basic quinoline pattern
    quinoline_matches = mol.GetSubstructMatches(quinoline_pattern)
    
    if not quinoline_matches:
        # Try alternative quinoline pattern (can catch some variations)
        quinoline_pattern2 = Chem.MolFromSmarts('c1cccc2c1[n]c(cc2)')
        quinoline_matches = mol.GetSubstructMatches(quinoline_pattern2)
        if not quinoline_matches:
            return False, "No quinoline core structure found"

    # Check for aromaticity
    aromatic_atoms = [atom.GetIsAromatic() for atom in mol.GetAtoms()]
    if not any(aromatic_atoms):
        return False, "Molecule is not aromatic"

    return True, "Contains quinoline core and one carboxylic acid group"
# Pr=None
# Recall=None