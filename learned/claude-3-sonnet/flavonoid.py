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

    # Basic benzopyran core pattern (O1c2ccccc2CC1)
    benzopyran_pattern = Chem.MolFromSmarts("O1c2ccccc2CC1")
    
    # Chromene core pattern (O1c2ccccc2C=C1)
    chromene_pattern = Chem.MolFromSmarts("O1c2ccccc2C=C1")
    
    # Chromenone (4-oxo) pattern
    chromenone_pattern = Chem.MolFromSmarts("O1c2ccccc2C(=O)C=C1")
    
    # Check for any of the core patterns
    has_benzopyran = mol.HasSubstructMatch(benzopyran_pattern)
    has_chromene = mol.HasSubstructMatch(chromene_pattern)
    has_chromenone = mol.HasSubstructMatch(chromenone_pattern)
    
    if not (has_benzopyran or has_chromene or has_chromenone):
        return False, "No benzopyran/chromene core structure found"

    # Look for aromatic ring attached at position 2
    # Multiple SMARTS to catch different variations
    aryl_patterns = [
        Chem.MolFromSmarts("O1[cX3]2[cX3][cX3][cX3][cX3][cX3]2C(Ar)=C1"),  # chromene with aryl
        Chem.MolFromSmarts("O1[cX3]2[cX3][cX3][cX3][cX3][cX3]2C(Ar)C1"),   # benzopyran with aryl
        Chem.MolFromSmarts("O1[cX3]2[cX3][cX3][cX3][cX3][cX3]2C(Ar)C(=O)C1")  # flavanone pattern
    ]
    
    has_aryl = any(mol.HasSubstructMatch(pattern) for pattern in aryl_patterns if pattern is not None)
    if not has_aryl:
        return False, "No aryl substituent at position 2"

    # Common substitution patterns in flavonoids
    # Look for hydroxyl or methoxy groups
    oh_pattern = Chem.MolFromSmarts("[OH]")
    ome_pattern = Chem.MolFromSmarts("[OX2][CH3]")
    
    num_oh = len(mol.GetSubstructMatches(oh_pattern)) if oh_pattern else 0
    num_ome = len(mol.GetSubstructMatches(ome_pattern)) if ome_pattern else 0
    
    # Most flavonoids have at least one OH or OMe group
    if num_oh + num_ome == 0:
        return False, "No typical flavonoid substitution pattern found"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Additional check for typical atom counts
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:  # Flavonoids typically have at least 15 atoms
        return False, "Molecule too small to be a flavonoid"

    return True, "Contains benzopyran core with aryl substituent at position 2 and typical flavonoid substitution pattern"