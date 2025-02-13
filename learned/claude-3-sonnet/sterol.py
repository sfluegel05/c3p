"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:36662 sterol

A sterol is defined as any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
(additional carbon atoms may be present in the side chain).

Examples:
- cholesterol
- ergosterol
- stigmasterol
- campesterol

"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = len(mol.GetSubstructMatches(hydroxy_pattern))
    if hydroxy_matches != 1:
        return False, f"Found {hydroxy_matches} hydroxy groups, need exactly 1"
    
    # Check for steroid backbone
    steroid_patterns = [
        Chem.MolFromSmarts("[C@]12[C@H]([C@@H](CC1)[C@@H](CC2)[H])[H]"), # cholestane
        Chem.MolFromSmarts("[C@@]12[C@H](CC[C@H]1[C@@H](CC2)[H])[H]"),   # cholane
        Chem.MolFromSmarts("[C@@]12[C@H]([C@@H](CC1)[C@@H](CC2)[H])[H]"), # ergostane
        Chem.MolFromSmarts("[C@@]12[C@H](CC[C@H]1[C@@H](CC2)[H])[H]")     # ergostane
    ]
    
    backbone_match = any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns)
    if not backbone_match:
        return False, "No steroid backbone found"
    
    # Check for cyclopropane ring
    cyclopropane_pattern = Chem.MolFromSmarts("[C@@]1([C@]2([C@@]1([H])[H])[C@@]2([H])[H])[H]")
    cyclopropane_match = mol.HasSubstructMatch(cyclopropane_pattern)
    
    # Check for side chains
    side_chain_pattern = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    
    # Combine criteria
    if backbone_match and (cyclopropane_match or len(side_chain_matches) > 0):
        return True, "Contains steroid backbone and meets other criteria for sterols"
    else:
        return False, "Does not meet structural criteria for sterols"