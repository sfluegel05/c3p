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
    # [o+] containing heterocycle fused to an aromatic ring and connected to another aromatic ring
    flavylium_core = Chem.MolFromSmarts('[o+]1c(-[c])c([cH,cO,cC])c([cH,cO,cC])c2c([cH,cO,cC])c([cH,cO,cC])c([cH,cO,cC])cc12')
    if not mol.HasSubstructMatch(flavylium_core):
        return False, "No flavylium cation core structure found"

    # Must have at least two hydroxyl groups on the main ring system
    hydroxyl_pattern = Chem.MolFromSmarts('(c(O)cc(O))||(c(O)c(O))')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic hydroxyl pattern"

    # Count oxygen atoms (excluding the charged oxygen)
    oxygen_count = sum(1 for atom in mol.GetAtoms() 
                      if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0)
    if oxygen_count < 2:
        return False, "Insufficient oxygen substituents"

    # Check for phenyl ring attachment
    phenyl_pattern = Chem.MolFromSmarts('c1ccccc1')
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "Missing phenyl ring substituent"

    # Common substituent patterns in anthocyanidins
    substituents = {
        'hydroxyl': Chem.MolFromSmarts('cO[H]'),
        'methoxy': Chem.MolFromSmarts('cOC'),
        'glycoside': Chem.MolFromSmarts('OC1OC(CO)C(O)C(O)C1O'),
        'acyl': Chem.MolFromSmarts('C(=O)'),
    }
    
    found_substituents = []
    for name, pattern in substituents.items():
        if mol.HasSubstructMatch(pattern):
            found_substituents.append(name)

    # Build classification message
    base_message = "Contains flavylium cation core"
    if found_substituents:
        substituents_str = ", ".join(found_substituents)
        message = f"{base_message} with {substituents_str} substituents"
    else:
        message = base_message

    return True, message