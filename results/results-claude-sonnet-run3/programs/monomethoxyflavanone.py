from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monomethoxyflavanone(smiles: str):
    """
    Determines if a molecule is a monomethoxyflavanone - a flavanone substituted by a single methoxy group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monomethoxyflavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for flavanone core structure (chroman-4-one)
    flavanone_pattern = Chem.MolFromSmarts('O1c2ccccc2C(=O)CC1')
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Not a flavanone core structure"
        
    # Check for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts('OC')
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    # Count methoxy groups
    methoxy_count = 0
    for match in methoxy_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        
        # Check if oxygen is connected to aromatic carbon
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                methoxy_count += 1
                break
    
    if methoxy_count == 1:
        return True, "Contains flavanone core with single methoxy substituent"
    elif methoxy_count == 0:
        return False, "No methoxy groups found"
    else:
        return False, f"Found {methoxy_count} methoxy groups (only one allowed)"
# Pr=None
# Recall=0.0