"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium (2-phenylchromenylium).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for positive charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 1:
        return False, "Must have a total charge of +1"
        
    # Look for flavylium core structure (2-phenylchromenylium)
    # [O+]=C1C=CC=CC=C1 connected to phenyl ring
    flavylium_pattern = Chem.MolFromSmarts('[O+]=C1c2ccccc2C=Cc3ccccc13')
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium core structure found"
        
    # Must have oxygen atoms (typically hydroxyls) on the rings
    oxygen_pattern = Chem.MolFromSmarts('cO')
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 2:
        return False, "Insufficient hydroxyl/oxygen substituents"
        
    # Check for characteristic benzopyran structure
    benzopyran_pattern = Chem.MolFromSmarts('c1c2c(cc1)OC=CC2')
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran structure found"
        
    # Count number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Must have at least 3 rings"
        
    # Look for characteristic hydroxyl pattern on A and B rings
    hydroxyl_pattern = Chem.MolFromSmarts('c(O)c(O)')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic hydroxyl pattern"
        
    # Additional checks for substituents commonly found in anthocyanidins
    # (optional methoxy groups)
    methoxy_pattern = Chem.MolFromSmarts('cOC')
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    
    # Success message varies based on presence of methoxy groups
    base_message = "Contains flavylium core with characteristic hydroxyl pattern"
    if has_methoxy:
        message = base_message + " and methoxy substituents"
    else:
        message = base_message
        
    return True, message