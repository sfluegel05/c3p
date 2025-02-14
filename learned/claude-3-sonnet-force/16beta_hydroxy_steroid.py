"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:132026 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is a steroid with a hydroxy group at position 16 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]1(CC[C@]2([C@@]1([C@H]([C@@]3([C@H]([C@@H]2[C@H](C3)C)C)C)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 16beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([C@@H]1[C@H]2[C@@H]([C@@H]([C@H](C2)[C@@H](O)[C@H]1C)C)C)O")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 16beta-hydroxy group found"
    
    # Check for common steroid ring systems and double bond patterns
    ring_systems = ["[C@]12[C@@H]([C@@]3([C@H]([C@@H]1[C@H](C2)C)CC[C@@H]4[C@@]3(CCC[C@H](C4)O)C)C)",
                    "[C@]12[C@@H]([C@@]3([C@H]([C@@H]1[C@H](C2)C)CC[C@@H]4[C@@]3(CC[C@H](C4)O)C)C)"]
    has_ring_system = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in ring_systems)
    if not has_ring_system:
        return False, "Steroid ring system not recognized"
    
    # Check for common double bond patterns
    double_bond_patterns = ["[C@@H]1[C@H](C=C)[C@@H]2[C@@H]3[C@H](C=C)[C@H](C=C)[C@H](C2)C31C",
                            "[C@@H]1[C@H](C=C)[C@@H]2[C@@H]3[C@H](C=C)[C@H](C=C)[C@H](C2)[C@H]31"]
    has_double_bonds = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in double_bond_patterns)
    if not has_double_bonds:
        return False, "Steroid double bond pattern not recognized"
    
    return True, "Contains steroid backbone with 16beta-hydroxy group"