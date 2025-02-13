"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36506 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern
    steroid_backbone_pattern = Chem.MolFromSmarts("[C@]12CCC[C@]3([C@@H]4CC[C@@]([C@](C)(CC[C@]2([H])CC1)(C)C4)C)C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
    
    # Look for 3beta-hydroxyl group pattern
    hydroxy_pattern = Chem.MolFromSmarts("[C@@]1(CC[C@@]2([H])CC[C@]3([H])C[C@@H](O)CC[C@]23C)[C@H](C)CCCC")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxyl group found"
    
    # Check for additional structural features
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for steroid"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms for steroid"
    
    # Count rings to ensure tetracyclic structure
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings != 4:
        return False, "Not a tetracyclic structure"
    
    return True, "Contains steroid backbone and 3beta-hydroxyl group"