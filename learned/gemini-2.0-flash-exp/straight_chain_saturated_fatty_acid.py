"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Reason for the classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Verify Carboxylic Acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
         return False, "No carboxylic acid group found"
    
    # Check for saturation
    if not all(bond.GetBondType() == Chem.BondType.SINGLE for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6):
        return False, "Not a saturated molecule"

    # Check for straight chain (no branching)
    branch_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])")
    if mol.HasSubstructMatch(branch_pattern):
      return False, "Branched chain detected"

    # Check for long carbon chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbon atoms"
      
    # Check for zero or one hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    n_hydroxy = len(mol.GetSubstructMatches(hydroxy_pattern))
    if n_hydroxy > 1:
        return False, "Too many hydroxy groups"
    
    # Check number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < (c_count - 2) or n_rotatable > c_count: # to deal with the optional OH group
        return False, f"Wrong number of rotatable bonds {n_rotatable}, number of C atoms {c_count}"

    return True, "Straight-chain saturated fatty acid"