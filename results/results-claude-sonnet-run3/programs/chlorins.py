from rdkit import Chem
from rdkit.Chem import AllChem

def is_chlorins(smiles: str):
    """
    Determines if a molecule is a chlorin (dihydroporphyrin with saturated carbons 
    in one pyrrole ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorin, False otherwise
        str: Reason for classification
    """
    # Handle molecules with Mg ion by removing it
    if '.' in smiles:
        smiles = smiles.split('.')[-1]
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for the core porphyrin-like structure using SMARTS patterns
    # Pattern for pyrrole rings (both aromatic and non-aromatic)
    pyrrole_pattern = Chem.MolFromSmarts('[nH,n-]1[cH,c][cH,c][cH,c][cH,c]1')
    
    # Pattern for reduced pyrrole ring (with two saturated carbons)
    reduced_pyrrole = Chem.MolFromSmarts('[nH,n-]1[CH2][CH2][cH,c][cH,c]1')
    
    # Check for presence of pyrrole rings
    matches = mol.GetSubstructMatches(pyrrole_pattern)
    reduced_matches = mol.GetSubstructMatches(reduced_pyrrole)
    
    total_rings = len(matches) + len(reduced_matches)
    
    if total_rings < 4:
        return False, "Does not contain enough pyrrole rings for chlorin structure"
        
    if len(reduced_matches) == 0:
        return False, "No reduced pyrrole ring found"
        
    if len(reduced_matches) > 1:
        return False, "More than one reduced pyrrole ring found"

    # Pattern to check for macrocyclic connectivity
    macrocycle_pattern = Chem.MolFromSmarts('[nH,n-]1[cH,c][cH,c][cH,c][cH,c]1-[cH,c]2[cH,c][cH,c][cH,c][nH,n-]2')
    
    if not mol.HasSubstructMatch(macrocycle_pattern):
        return False, "Pyrrole rings not connected in macrocyclic structure"

    # Check for conjugated system
    conjugated_system = Chem.MolFromSmarts('c1cccc(c1)-c2ccccc2')
    if not mol.HasSubstructMatch(conjugated_system):
        return False, "Missing required conjugated system"

    return True, "Molecule contains chlorin structure with one reduced pyrrole ring"
# Pr=None
# Recall=0.0