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
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@@H]([C@@]3([C@H]([C@@H]1[C@H](C2)C)CC[C@@H]4[C@@]3(CC[C@H](C4)O)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 16beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([C@H]1[C@@H]2[C@H]([C@@H]([C@@H]([C@H](C2)C)O)[C@H]1C)C)O")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 16beta-hydroxy group found"
    
    # Check molecular weight - steroids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for steroid"
    
    return True, "Contains steroid backbone with 16beta-hydroxy group"