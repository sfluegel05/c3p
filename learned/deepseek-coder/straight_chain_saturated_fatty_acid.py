"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: CHEBI:26607 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is a carboxylic acid with a straight carbon chain and no double/triple bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for unsaturation (double/triple bonds) in the carbon chain
    # Exclude the carboxylic acid group from the unsaturation check
    for bond in mol.GetBonds():
        if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.AROMATIC]:
            # Check if the bond is part of the carboxylic acid group
            if not (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE):
                return False, "Molecule contains double or triple bonds (unsaturated)"

    # Check for branching in the carbon chain
    # We look for carbons connected to more than 2 other carbons, excluding the carboxylic acid group
    carbon_chain_branch_pattern = Chem.MolFromSmarts("[CX4;H0,H1,H2]([CX4])([CX4])[CX4]")
    if mol.HasSubstructMatch(carbon_chain_branch_pattern):
        return False, "Molecule contains branches in the carbon chain"

    # Check carbon chain length (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short (minimum 4 carbons required)"

    return True, "Straight-chain saturated fatty acid with a carboxylic acid group"