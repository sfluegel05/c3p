"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
A carbohydrate derivative where one or more oxygens/hydroxy groups are replaced by sulfur or -SR
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Must contain sulfur
    if not any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "No sulfur atoms present"

    # Look for pyranose or furanose ring patterns
    sugar_patterns = [
        "[C]1[O,S][C]([C,O,S])[C]([O,S])[C]([O,S])[C]1[O,S]",  # pyranose
        "[C]1[O,S][C]([C,O,S])[C]([O,S])[C]1[O,S]"  # furanose
    ]
    
    found_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_sugar = True
            break
            
    if not found_sugar:
        return False, "No sugar ring structure found"

    # Count carbons, oxygens and sulfurs in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # Basic sanity checks for sugar-like composition
    if c_count < 4:
        return False, "Too few carbons for a sugar structure"
    
    if o_count + s_count < 4:
        return False, "Too few O/S atoms for a sugar structure"

    # Look specifically for sulfur connected to the sugar ring
    sugar_s_pattern = "[C]1[O][C]([!O])[C]([O,S])[C]([O,S])[C]1[O,S]"
    sugar_with_s = mol.HasSubstructMatch(Chem.MolFromSmarts(sugar_s_pattern))
    
    # Additional check for S-glycosidic bond
    s_glycosidic = mol.HasSubstructMatch(Chem.MolFromSmarts("[C]1[O][C]([S])[C]([O])[C]([O])[C]1[O]"))
    
    if not (sugar_with_s or s_glycosidic):
        return False, "No sulfur atoms connected to sugar structure"

    # Check for characteristic sugar hydroxyl/thiol pattern
    hydroxy_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[C][O,S][H]")))
    if hydroxy_count < 2:
        return False, "Too few hydroxyl/thiol groups for a sugar derivative"

    return True, "Contains sugar ring with sulfur substitution"