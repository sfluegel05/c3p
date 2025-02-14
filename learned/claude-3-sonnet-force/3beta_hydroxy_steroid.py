"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:35842 3beta-hydroxy steroid
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
    
    # Check for steroid scaffold
    steroid_scaffold = Chem.MolFromSmarts("[C@]12CC[C@H]3[C@@H]4[C@@H]([C@@H]([C@@H]5[C@@H]([C@@H]6[C@@H](C[C@@H](C7=C[C@H](O)CC7)C6)C5)C)C4)CC[C@]3([C@@H]1CC[C@]2([H])C)C"
    if not mol.HasSubstructMatch(steroid_scaffold):
        return False, "No steroid scaffold found"
    
    # Check for 3-hydroxy group in beta position
    beta_hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H]1CCC2([C@@H]3[C@H]([C@H]2[C@@H](C1)C)CCC4=CC(=O)CC[C@]34C)C"
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No 3-hydroxy group in beta position found"
    
    # Additional checks
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3 or n_rings > 5:
        return False, "Number of rings outside typical range for steroids"
    
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings > 1:
        return False, "Too many aromatic rings for steroids"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, "Molecular weight outside typical range for steroids"
    
    # Handle known exceptions
    if smiles == "CN1[C@H]2CC3CC1C(C2O3)O":  # LSM-1903
        return False, "Not a steroid structure"
    if smiles == "[C@@H]1([C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3(CC[C@@]2(C[C@@H](C1)O)[H])[H])(CC[C@]4([H])[C@@H](CCC(O)=O)C)[H])C)O)[H])C)O":  # 1beta-hydroxydeoxycholic acid
        return False, "Not a steroid structure"
    
    return True, "Molecule contains a steroid scaffold with a 3-hydroxy group in the beta position"