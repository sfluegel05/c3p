"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is characterized by the cyclopenta[a]phenanthrene carbon skeleton,
    with methyl groups at C-10 and C-13, and often an alkyl group at C-17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid core pattern (cyclopenta[a]phenanthrene skeleton)
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C@@]2([C@@]([C@@]3([C@@]([C@@]4([C@](CC3)(C[C@H](CC4)C)C)[H])[H])(CC2)[H])[H])CC1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No cyclopenta[a]phenanthrene skeleton found"

    # Check for methyl groups at C-10 and C-13
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 2:
        return False, "Missing methyl groups at C-10 and/or C-13"

    # Check for an alkyl group at C-17 (optional, but common)
    alkyl_pattern = Chem.MolFromSmarts("[CX4]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    if len(alkyl_matches) < 3:
        return False, "Missing alkyl group at C-17"

    # Count carbons and hydrogens to ensure a reasonable steroid structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if c_count < 17:
        return False, "Too few carbons for a steroid"
    if h_count < 28:
        return False, "Too few hydrogens for a steroid"

    return True, "Contains cyclopenta[a]phenanthrene skeleton with characteristic methyl and alkyl groups"