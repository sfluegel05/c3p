"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: CHEBI:38043 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is a steroid with a hydroxy group at position 11 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]12CCC[C@H]3[C@@H]4CC[C@@]5(C)[C@@H](O)CC[C@]5(O)[C@H]4CC[C@]3(C)C1=O")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 11beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C@@]1")
    hydroxy_match = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_match:
        return False, "No 11beta-hydroxy group found"
    
    # Check for additional features
    n_ketones = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()) == 1)
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    if n_ketones < 1 or n_double_bonds < 1 or mol_wt < 250 or mol_wt > 600:
        return False, "Missing essential features for 11beta-hydroxy steroid"
    
    return True, "Contains steroid backbone with 11beta-hydroxy group and essential features"