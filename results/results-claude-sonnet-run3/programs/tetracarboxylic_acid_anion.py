from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetracarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a tetracarboxylic acid anion.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tetracarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Count carboxylate groups (-C([O-])=O)
    carboxylate_pattern = Chem.MolFromSmarts('[C](=[O])[O-]')
    carboxylate_matches = len(mol.GetSubstructMatches(carboxylate_pattern))
    
    # Count carboxylic acid groups (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts('[C](=[O])[OH]')
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    total_carboxy = carboxylate_matches + carboxyl_matches
    
    if total_carboxy < 4:
        return False, f"Contains only {total_carboxy} carboxy groups (minimum 4 required)"
        
    if carboxylate_matches == 0:
        return False, "No deprotonated carboxy groups found"
        
    if total_carboxy == 4:
        if carboxylate_matches == 4:
            return True, "Fully deprotonated tetracarboxylic acid anion"
        else:
            return True, f"Partially deprotonated tetracarboxylic acid anion ({carboxylate_matches} of 4 groups deprotonated)"
            
    return True, f"Polycarboxylic acid anion with {total_carboxy} carboxy groups, {carboxylate_matches} deprotonated"
# Pr=None
# Recall=None