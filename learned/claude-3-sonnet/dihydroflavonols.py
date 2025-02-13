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

    # Basic dihydroflavonol core:
    # - Chromanone core (benzopyran-4-one)
    # - OH at position 3
    # - Phenyl at position 2
    # Note: Using more flexible pattern that allows substitutions
    core_pattern = Chem.MolFromSmarts('[O;!H0]-[CH1]-1-[c]2[c]([c,n])[c,n][c,n][c,n]2-C(=O)-[CH]([O;!H0])-[CH1]-1-[c]1[c][c][c][c][c]1')
    
    if core_pattern is None:
        return None, "Invalid SMARTS pattern"
    
    if not mol.HasSubstructMatch(core_pattern):
        # Try alternative pattern allowing for any substitution at position 2
        alt_pattern = Chem.MolFromSmarts('[O]1[CH1][c]2[c]([c,n])[c,n][c,n][c,n]2C(=O)[CH]([O;!H0])[CH1]1[c]1[c,n][c,n][c,n][c,n][c,n]1')
        if alt_pattern is None:
            return None, "Invalid alternative SMARTS pattern"
        if not mol.HasSubstructMatch(alt_pattern):
            return False, "Missing dihydroflavonol core structure"

    # Count basic features to ensure reasonable composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for dihydroflavonol structure"
    if o_count < 3:
        return False, "Too few oxygens for dihydroflavonol structure"

    # Check for ring systems
    ri = mol.GetRingInfo()
    if ri.NumRings() < 3:
        return False, "Missing required ring systems"

    return True, "Contains dihydroflavonol core structure with 3-OH group"