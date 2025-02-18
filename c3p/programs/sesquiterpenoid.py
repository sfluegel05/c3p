"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a 15-carbon skeleton derived from three isoprene units,
    and can have modifications such as removal of methyl groups, or additional functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Carbon count - must have at least 15 carbons, up to 30 to account for functional groups.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (15 <= c_count <= 30):
        return False, f"Incorrect number of carbons: {c_count}, must be between 15 and 30."

    # 2. Check for isoprene units (C5H8)
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C[C]") # Check for isoprene
    isopentane_pattern = Chem.MolFromSmarts("CC(C)CC") # Check for isopentane fragments, more common in modified molecules
    
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    isopentane_matches = mol.GetSubstructMatches(isopentane_pattern)
    
    if len(isoprene_matches) + len(isopentane_matches) < 2 : # require at least 2 isoprene or isopentane fragments
        return False, f"Too few isoprene/isopentane units: {len(isoprene_matches) + len(isopentane_matches)}, need at least 2"
    
    # 3. Check for common sesquiterpenoid frameworks - focus on rings with a size of 5,6,7 or open chains.
    # This is more flexible than a single very specific SMARTS.
    
    core_patterns = [
        Chem.MolFromSmarts("[C;R5]1[C;R5][C;R5][C;R5][C;R5]1"),  # 5-membered ring
        Chem.MolFromSmarts("[C;R6]1[C;R6][C;R6][C;R6][C;R6][C;R6]1"), # 6-membered ring
        Chem.MolFromSmarts("[C;R7]1[C;R7][C;R7][C;R7][C;R7][C;R7][C;R7]1"), #7-membered ring
        Chem.MolFromSmarts("[C;R0]([C;R0])([C;R0])([C;R0])[C;R0][C;R0]([C;R0])([C;R0])([C;R0])[C;R0][C;R0]([C;R0])([C;R0])[C;R0]"), #15 carbon chain
        Chem.MolFromSmarts("[C;R1]1[C;R1]2[C;R1]3[C;R1]1[C;R1]2[C;R1]3"), # a tricyclic structure
    ]
    
    found_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            found_core = True
            break
    
    if not found_core:
         return False, "No common sesquiterpenoid skeleton found"


    return True, "Contains a sesquiterpenoid core with at least 2 isoprene units."