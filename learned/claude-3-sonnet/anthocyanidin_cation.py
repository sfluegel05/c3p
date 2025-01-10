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

    # Look for the basic flavylium cation core:
    # More flexible pattern for the chromenylium system with phenyl substituent
    flavylium_core = Chem.MolFromSmarts('[o+]1c(c2)cc(c1)-c3ccccc3')
    if not mol.HasSubstructMatch(flavylium_core):
        return False, "No flavylium cation core structure found"

    # Check for the fused benzene ring
    fused_ring = Chem.MolFromSmarts('[o+]1c(c2ccccc2)cc*c1')
    if not mol.HasSubstructMatch(fused_ring):
        return False, "Missing required fused ring system"

    # Count oxygen-containing substituents (excluding the charged oxygen)
    oxygen_pattern = Chem.MolFromSmarts('cO')
    oxygen_matches = len(mol.GetSubstructMatches(oxygen_pattern))
    if oxygen_matches < 2:
        return False, "Insufficient oxygen substituents"

    # Look for characteristic substitution patterns
    patterns = {
        'hydroxyl': Chem.MolFromSmarts('cO[H]'),
        'methoxy': Chem.MolFromSmarts('cOC'),
        'glycoside': Chem.MolFromSmarts('OC1OC(CO)C(O)C'),
        'acyl': Chem.MolFromSmarts('C(=O)'),
    }
    
    found_patterns = []
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            found_patterns.append(name)
            
    if not found_patterns:
        return False, "Missing characteristic oxygen-containing substituents"

    # Additional check for proper ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring systems"

    # Build classification message
    base_msg = "Contains flavylium cation core with proper ring system"
    if found_patterns:
        substituents_str = ", ".join(found_patterns)
        message = f"{base_msg} and {substituents_str} substituents"
    else:
        message = base_msg

    return True, message