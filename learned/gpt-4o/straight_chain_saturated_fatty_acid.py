"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

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

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Ensure all bonds in the main chain are single bonds (saturated)
    if not all(bond.GetBondType() == Chem.BondType.SINGLE for bond in mol.GetBonds()):
        return False, "Presence of double or triple bonds indicates unsaturation"

    # Determine if molecule is a straight chain
    # Count carbon atoms in the longest chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # We assume that a fatty acid will have n+1 carbons when saturated with a terminal carboxylic group
    if mol.GetNumHeavyAtoms() != c_count + 2:
        return False, "Presence of branch or other functional group"

    # Hydroxy group allowance: count OH groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))
    if hydroxy_count > 1:
        return False, f"Too many hydroxy groups: found {hydroxy_count}, need at most 1"

    return True, "Molecule is a straight-chain saturated fatty acid"

# Example use:
# smiles = "CCCCCCCC(O)=O"  # Octanoic acid
# result, reason = is_straight_chain_saturated_fatty_acid(smiles)
# print(result, reason)