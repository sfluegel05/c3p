"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: CHEBI:26815 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is a steroid with a hydroxyl group in the beta configuration
    at position 3, and a double bond between positions 5 and 6 in the Delta(5) configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_backbone_pattern = Chem.MolFromSmarts("[C@@]12[C@H]([C@@]([C@H](C1)C)([C@@]2([H])[H])[H])[H]")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
    
    # Check for hydroxyl group at position 3 in beta configuration
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@]([C@H](O)[H])(C)"
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No 3-beta hydroxyl group found"
    
    # Check for double bond between positions 5 and 6 in Delta(5) configuration
    delta5_pattern = Chem.MolFromSmarts("[C@@]1([H])=C([C@@H](C)C)[C@@H](C)[C@H](C)C1")
    if not mol.HasSubstructMatch(delta5_pattern):
        return False, "No Delta(5) double bond found"
    
    # Additional checks
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Too few rotatable bonds for a steroid"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Too few carbons for a steroid"
    
    # If all checks pass, classify as a 3beta-hydroxy-Delta(5)-steroid
    return True, "Contains a steroid backbone with a 3-beta hydroxyl group and a Delta(5) double bond"