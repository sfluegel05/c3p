"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is a fatty acid with a hydroxyl group on the carbon
    in the beta- or 3-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a hydroxy group at the 3rd carbon position from the carboxylic acid group
    beta_hydroxy_pattern = Chem.MolFromSmarts("C(=O)[O]-[C]-[C]-[C]-[OH]")
    hydroxy_matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
    if not hydroxy_matches:
         return False, "No hydroxy group found at beta position"

    #Check that the hydroxyl is part of a long chain containing at least one additional carbon past the hydroxyl
    long_chain_pattern = Chem.MolFromSmarts("[C]([OH])-[C]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if not long_chain_matches:
        return False, "Hydroxy group is not on a long chain"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 4 or o_count < 3:
        return False, "Too few carbons or oxygens for fatty acid"

    # Count rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 1:
         return False, "Too few rotatable bonds for fatty acid"

    return True, "3-Hydroxy fatty acid identified"