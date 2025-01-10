"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are characterized by a steroid nucleus with a hydroxyl group 
    at C3 and possible additional alkyl groups in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C"),  # basic steroid skeleton
        Chem.MolFromSmarts("C1=CC[C@H]2[C@@H]3CC[C@@]4([H])CCC[C@]4(C)[C@H]3CC[C@@]21[H]")  # alternate stereo/unsaturated
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No steroid backbone found"
        
    # Check for hydroxyl group at position 3
    hydroxy_position_3 = Chem.MolFromSmarts("O[C@@H]1CC[C@H](C2CCCC3C2CCC4CC(C)(C)C3(C)CC4)C1")
    if not mol.HasSubstructMatch(hydroxy_position_3):
        return False, "Missing hydroxyl group at position 3"

    # Check for typical phytosterol double bond positions, e.g., at C5, C6
    double_bond_patterns = [
        Chem.MolFromSmarts("C1=CC[C@H]2[C@@H]3CC=C4C(C2)[C@H]3CC[C@@]41C")  # common phytosterol unsaturations
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in double_bond_patterns):
        return False, "No expected double bonds found in rings"

    # Check for side chain modifications (e.g., methyl or ethyl groups)
    side_chain_patterns = [
        Chem.MolFromSmarts("[CX4,CX3](C)"),  # broader search for alkylation
        Chem.MolFromSmarts("[C@H](C)[C@H](CC)C(C)C")  # specific for C24 modification linked to sitosterol-like structures
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in side_chain_patterns):
        return False, "No typical phytosterol side chain alkylations found"

    return True, "Steroid nucleus with typical phytosterol hydroxyl, double bonds, and side chain modifications found"