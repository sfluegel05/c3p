from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinolinemonocarboxylate(smiles: str):
    """
    Determines if a molecule is a quinoline monocarboxylate.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinoline monocarboxylate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for presence of carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts('[C](=[O])[O-]')
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if not carboxylate_matches:
        return False, "No carboxylate group found"
    
    if len(carboxylate_matches) > 1:
        return False, "More than one carboxylate group found"
        
    # Check for quinoline core structure
    quinoline_pattern = Chem.MolFromSmarts('c1cccc2c1[n]ccc2')
    quinoline_matches = mol.GetSubstructMatches(quinoline_pattern)
    
    if not quinoline_matches:
        # Try alternative quinoline pattern with different atom ordering
        quinoline_pattern2 = Chem.MolFromSmarts('c1ccc2c(c1)ccc[n]2')
        quinoline_matches = mol.GetSubstructMatches(quinoline_pattern2)
        
        if not quinoline_matches:
            return False, "No quinoline core structure found"
    
    # Verify carboxylate is attached to quinoline
    carboxylate_carbon = carboxylate_matches[0][0]
    quinoline_atoms = set(sum(quinoline_matches, ()))
    
    # Get neighbors of carboxylate carbon
    neighbors = [x.GetIdx() for x in mol.GetAtomWithIdx(carboxylate_carbon).GetNeighbors()]
    
    # Check if any neighbor is part of quinoline core
    if not any(n in quinoline_atoms for n in neighbors):
        return False, "Carboxylate group not attached to quinoline core"
        
    return True, "Valid quinoline monocarboxylate structure found"
# Pr=1.0
# Recall=0.8333333333333334