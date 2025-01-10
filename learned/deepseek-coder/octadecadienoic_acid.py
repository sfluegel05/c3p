"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:32395 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain, C18 polyunsaturated fatty acid with two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Expected 18 carbon atoms, found {c_count}"

    # Count the number of double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count != 2:
        return False, f"Expected 2 double bonds, found {double_bond_count}"

    # Check if the molecule is a straight-chain fatty acid
    # This is a heuristic check: we expect the molecule to have a long carbon chain
    # with no branching (i.e., no more than 2 bonds per carbon except for the ends)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if atom.GetDegree() > 2 and not (atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3):  # End of chain
                return False, "Molecule is not a straight-chain fatty acid"

    return True, "Straight-chain C18 fatty acid with 2 double bonds and a carboxylic acid group"