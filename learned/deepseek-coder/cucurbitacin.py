"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: CHEBI:39448 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for tetracyclic structure (4 rings)
    n_rings = len(Chem.GetSSSR(mol))
    if n_rings < 4:
        return False, f"Found {n_rings} rings, need at least 4 for tetracyclic structure"

    # Check for triterpenoid skeleton (30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Found {c_count} carbons, need at least 30 for triterpenoid skeleton"

    # Check for key functional groups (hydroxyl and carbonyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found"
    if len(carbonyl_matches) < 1:
        return False, "No carbonyl groups found"

    # Check for cucurbitane skeleton pattern (specific ring arrangement)
    cucurbitane_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]3[C][C][C]4[C]5[C][C][C]1[C]2[C]3[C]4[C]5")
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "No cucurbitane skeleton found"

    return True, "Contains tetracyclic triterpenoid structure with cucurbitane skeleton and key functional groups"