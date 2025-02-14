"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: CHEBI:35857 Methyl-branched fatty acid
A methyl-branched fatty acid is defined as any branched-chain fatty acid containing methyl branches only.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for methyl branches
    methyl_pattern = Chem.MolFromSmarts("CC")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)

    # Check for linear carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCC")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)

    # If no methyl branches or no linear carbon chains, not a methyl-branched fatty acid
    if not methyl_matches or not carbon_chain_matches:
        return False, "No methyl branches or linear carbon chains found"

    # Iterate over atoms to check for other types of branches
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2:  # Carbon with more than 2 neighbors
            neighbors = [mol.GetAtomWithIdx(nei).GetAtomicNum() for nei in atom.GetNeighbors()]
            if any(nei != 6 and nei != 1 for nei in neighbors):  # Neighbor is not carbon or hydrogen
                return False, "Found non-methyl branches"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be a fatty acid"

    return True, "Molecule contains a carboxylic acid group, methyl branches, and linear carbon chains"