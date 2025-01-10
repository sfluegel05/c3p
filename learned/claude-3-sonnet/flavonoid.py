"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid has a 1-benzopyran core with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavonoid core patterns - more flexible versions
    core_patterns = [
        # Basic chromene/flavone core (O1c2ccccc2C(=O)C=C1)
        Chem.MolFromSmarts("O1c2c([c,C][c,C][c,C][c,C]2)C(=O)C=C1"),
        
        # Flavanone core (O1c2ccccc2C(=O)CC1)
        Chem.MolFromSmarts("O1c2c([c,C][c,C][c,C][c,C]2)C(=O)CC1"),
        
        # Isoflavone core
        Chem.MolFromSmarts("O1c2c([c,C][c,C][c,C][c,C]2)C(=C)C(=O)C1"),
        
        # More general benzopyran core
        Chem.MolFromSmarts("O1c2c([c,C][c,C][c,C][c,C]2)[C,c][C,c]1"),
        
        # Flavonol core
        Chem.MolFromSmarts("O1c2c([c,C][c,C][c,C][c,C]2)C(O)=C(C1=O)")
    ]
    
    has_core = False
    for pattern in core_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_core = True
            break
    
    if not has_core:
        return False, "No flavonoid core structure found"

    # Check for aromatic ring system connected to the core
    # More flexible patterns for the aryl substituent
    aryl_patterns = [
        # General aryl group at position 2
        Chem.MolFromSmarts("O1[c,C]2[c,C][c,C][c,C][c,C][c,C]2[C,c]([a;r6])[C,c]1"),
        
        # Alternative connection pattern
        Chem.MolFromSmarts("O1[c,C]2[c,C][c,C][c,C][c,C][c,C]2[C,c](=O)[C,c]([a;r6])")
    ]
    
    has_aryl = False
    for pattern in aryl_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_aryl = True
            break

    if not has_aryl:
        return False, "Missing required aryl substituent"

    # Check for typical flavonoid characteristics
    # Look for hydroxyl, methoxy groups, or glycosidic linkages
    substitution_patterns = [
        Chem.MolFromSmarts("[OH]"),  # hydroxyl
        Chem.MolFromSmarts("[OX2][CH3]"),  # methoxy
        Chem.MolFromSmarts("[OX2][CH]"),  # possible glycosidic linkage
        Chem.MolFromSmarts("O=C")  # carbonyl groups
    ]
    
    substitution_count = sum(
        len(mol.GetSubstructMatches(pattern)) 
        for pattern in substitution_patterns 
        if pattern is not None
    )
    
    if substitution_count < 2:  # Most flavonoids have multiple substituents
        return False, "Insufficient number of typical flavonoid substituents"

    # Ring count check - flavonoids typically have at least 3 rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Size check
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:
        return False, "Molecule too small to be a flavonoid"

    return True, "Contains flavonoid core structure with appropriate substituents"