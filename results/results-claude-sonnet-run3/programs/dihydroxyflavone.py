from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroxyflavone(smiles: str):
    """
    Determines if a molecule is a dihydroxyflavone (flavone with exactly two hydroxy substituents).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a dihydroxyflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Basic flavone core pattern (more specific than previous)
    flavone_pattern = Chem.MolFromSmarts('O=C1C=C(Oc2c1cccc2)-c1ccccc1')
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No flavone core structure found"

    # Find hydroxyl groups attached to aromatic carbons
    hydroxy_pattern = Chem.MolFromSmarts('c-[OH]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count unique aromatic carbons with OH groups
    unique_hydroxy_carbons = set()
    for match in hydroxy_matches:
        c_atom = match[0]  # Get the carbon atom index
        unique_hydroxy_carbons.add(c_atom)
    
    # Get all methoxy groups to ensure they're not counted
    methoxy_pattern = Chem.MolFromSmarts('cOC')
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    # Remove any carbons that have methoxy groups
    for match in methoxy_matches:
        if match[0] in unique_hydroxy_carbons:
            unique_hydroxy_carbons.remove(match[0])
    
    # Count the remaining unique hydroxyl groups
    hydroxy_count = len(unique_hydroxy_carbons)
    
    if hydroxy_count != 2:
        return False, f"Found {hydroxy_count} hydroxy groups, need exactly 2"
    
    # Get positions of OH groups for reporting
    positions = sorted(list(unique_hydroxy_carbons))
    
    # Additional check to ensure hydroxyls are on the flavone rings
    flavone_match = mol.GetSubstructMatch(flavone_pattern)
    flavone_atoms = set(flavone_match)
    
    # Verify OH groups are attached to the flavone core
    for pos in positions:
        if pos not in flavone_atoms:
            return False, "Hydroxy groups must be on flavone core structure"
            
    return True, f"Dihydroxyflavone with hydroxy groups at positions {positions}"
# Pr=None
# Recall=0.0