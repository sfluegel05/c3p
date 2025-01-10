"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid with one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find the carboxylic acid group (C(=O)[O,H])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find hydroxy groups (not directly part of carboxylic acid)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][O][H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy groups found"
    
    # Find carbon chain
    atom_counts = {atom.GetAtomicNum(): 0 for atom in mol.GetAtoms()}
    for atom in mol.GetAtoms():
        atom_counts[atom.GetAtomicNum()] += 1

    carbon_count = atom_counts.get(6, 0)
    # Ensure we have at least 6 carbons in non-aromatic chains
    chain_pattern = Chem.MolFromSmarts("[*]CCCC[*]")
    if not mol.HasSubstructMatch(chain_pattern) or carbon_count < 6:
        return False, "Carbon chain too short or not linear/branched"

    return True, "Contains one or more hydroxy groups and a carboxylic acid group, with sufficient and relevant carbon chain length"