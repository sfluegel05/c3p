"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: CHEBI:35577 17beta-hydroxy steroid
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is a steroid with a hydroxy group at position 17 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C&r1,r2,r3,r4]1CCC2C3=CC=C4[C@@]5(CC[C@H](C5)C4)C[C@H]3[C@@H](C2)C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Find the 17-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@]([H])(O)[C@@]1([H])CC[C@]2([H])C3=CC=C4[C@@]5(CC[C@H](C5)C4)C[C@H]3[C@@H](C2)C1")
    matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not matches:
        return False, "No 17-hydroxy group in the beta configuration found"

    # Additional checks
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 500:
        return False, "Molecular weight out of typical range for 17beta-hydroxy steroids"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 18 or o_count < 2:
        return False, "Incorrect number of carbon or oxygen atoms for 17beta-hydroxy steroids"

    return True, "Molecule has a 17beta-hydroxy group on the steroid backbone"