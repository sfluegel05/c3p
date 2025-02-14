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
    steroid_pattern = Chem.MolFromSmarts("[C@]1([C@H](CC[C@@]2([C@@]1([C@H]([C@@]3([C@H]([C@@H]2[C@H](C3)C)C)C)C)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 16beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([C@H]1[C@@H]2[C@H]([C@@H]([C@@H]([C@H](C2)C)O)[C@H]1C)C)O")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 16beta-hydroxy group found"
    
    # Check for common steroid ring systems
    ring_systems = ["[C@]12[C@@H]([C@@]3([C@H]([C@@H]1[C@H](C2)C)CC[C@@H]4[C@@]3(CC[C@H](C4)O)C)C)",
                    "[C@]12[C@@H]([C@@]3([C@H]([C@@H]1[C@H](C2)C)CC[C@@H]4[C@@]3(CCC[C@H](C4)O)C)C)"]
    has_ring_system = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in ring_systems)
    if not has_ring_system:
        return False, "Steroid ring system not recognized"
    
    # Check atom counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 17 or c_count > 30:
        return False, "Incorrect carbon count for steroid"
    if o_count < 1 or o_count > 4:
        return False, "Incorrect oxygen count for steroid"
    
    # Check molecular weight - steroids typically 200-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight out of typical steroid range"
    
    return True, "Contains steroid backbone with 16beta-hydroxy group"