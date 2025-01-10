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
    # A benzopyrylium ring system (O+ containing heterocycle) connected to a phenyl ring
    flavylium_core = Chem.MolFromSmarts('[O+]1=Cc2c(-c3ccc[c,n]c3)cc(O)cc2C=C1')
    if not mol.HasSubstructMatch(flavylium_core):
        # Try alternative pattern that might catch other variants
        alt_flavylium = Chem.MolFromSmarts('[O+]1=Cc2c(-[c,n]3[c,n]c[c,n][c,n][c,n]3)cc(O)cc2C=C1')
        if not mol.HasSubstructMatch(alt_flavylium):
            return False, "No flavylium cation core structure found"

    # Must have oxygen atoms (typically hydroxyls) on the rings
    # Look for hydroxyls on both the A and C rings
    hydroxyl_pattern = Chem.MolFromSmarts('c1c(O)cc(O)cc1')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic hydroxyl pattern on rings"

    # Count oxygen atoms (excluding the charged oxygen)
    oxygen_count = sum(1 for atom in mol.GetAtoms() 
                      if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0)
    if oxygen_count < 2:
        return False, "Insufficient oxygen substituents"

    # Common substituent patterns in anthocyanidins
    substituents = {
        'methoxy': Chem.MolFromSmarts('cOC'),
        'glycoside': Chem.MolFromSmarts('OC1OCC(O)C(O)C1O'),
        'acyl': Chem.MolFromSmarts('C(=O)C'),
    }
    
    found_substituents = []
    for name, pattern in substituents.items():
        if mol.HasSubstructMatch(pattern):
            found_substituents.append(name)

    # Build classification message
    base_message = "Contains flavylium cation core with hydroxyl groups"
    if found_substituents:
        substituents_str = ", ".join(found_substituents)
        message = f"{base_message} and {substituents_str} substituents"
    else:
        message = base_message

    return True, message