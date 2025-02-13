"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:15864 dihydroflavonol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxy group at position 3
    of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core dihydroflavonol structure:
    # - Chromanone core with OH at position 3
    # - More permissive pattern allowing substitutions
    core_pattern = Chem.MolFromSmarts('[O;!H0]-[CH1]-1-[c]2[c,n]([O]-1)[c,n][c,n][c,n][c,n]2-[C;!H0](=O)')
    if core_pattern is None:
        return None, "Invalid SMARTS pattern for core structure"
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroflavonol core structure"

    # Check for 2-phenyl group
    phenyl_pattern = Chem.MolFromSmarts('O1[CH1]([CH1](O)[C](=O)c2c1[c,n][c,n][c,n][c,n]2)[c]1[c,n][c,n][c,n][c,n][c,n]1')
    if phenyl_pattern is None:
        return None, "Invalid SMARTS pattern for 2-phenyl group"
        
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "Missing required 2-phenyl substituent"

    # Count basic features to ensure reasonable composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for dihydroflavonol structure"
    if o_count < 4:
        return False, "Too few oxygens for dihydroflavonol structure"

    # Check for ketone group at position 4
    ketone_pattern = Chem.MolFromSmarts('O1[CH1][CH1](O)[C](=O)c2c1[c,n][c,n][c,n][c,n]2')
    if ketone_pattern is None:
        return None, "Invalid SMARTS pattern for ketone group"
        
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group at position 4"

    # Verify 3-OH group specifically
    oh3_pattern = Chem.MolFromSmarts('O1[CH1][CH1](O)[C](=O)c2c1[c,n][c,n][c,n][c,n]2')
    if oh3_pattern is None:
        return None, "Invalid SMARTS pattern for 3-OH group"
        
    if not mol.HasSubstructMatch(oh3_pattern):
        return False, "Missing hydroxyl group at position 3"

    return True, "Contains dihydroflavonol core structure with 3-OH group and 2-phenyl substituent"