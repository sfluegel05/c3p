"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:24474 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic pattern for isoflavonoid core:
    # - Benzene ring fused to oxygen-containing heterocycle
    # - Aryl group at position 3
    # - Allows for variations in saturation and substitution
    core_pattern = Chem.MolFromSmarts('[c,C]1[c,C][c,C][c,C][c,C][c,C]2O[c,C][c,C]([$(c3ccccc3),$(C3=CC=CC=C3)])[c,C](=O)[c,C]12')
    
    # Alternative pattern for more complex fused systems
    alt_pattern = Chem.MolFromSmarts('[c,C]1[c,C][c,C][c,C][c,C]2O[c,C][c,C]([$(c3ccccc3),$(C3=CC=CC=C3)])[c,C](=O)[c,C]12')
    
    if core_pattern is None or alt_pattern is None:
        return None, "Error in SMARTS patterns"

    # Check for required substructures
    if not (mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(alt_pattern)):
        return False, "No isoflavonoid core structure found"

    # Verify presence of oxygen-containing heterocycle
    pyran_pattern = Chem.MolFromSmarts('c1c2c(cc1)OCC(=O)c2')
    if not mol.HasSubstructMatch(pyran_pattern):
        return False, "Missing benzopyran core"

    # Check for aryl group
    aryl_pattern = Chem.MolFromSmarts('a1aaaa1')
    aryl_matches = len(mol.GetSubstructMatches(aryl_pattern))
    if aryl_matches < 2:  # Need at least 2 aromatic rings (one from core, one from substituent)
        return False, "Missing required aryl substituent"

    # Count oxygens to ensure we have enough for the core structure
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Insufficient oxygen atoms for isoflavonoid structure"

    # Additional check for ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for isoflavonoid structure"

    return True, "Contains isoflavonoid core structure with aryl substituent at position 3"