"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid that contains three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")  # Include carboxylate ion as C(=O)[O-]
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # Count the number of C=C double bonds in an acyclic carbon chain
    # assuming linear or near-linear configuration
    double_bond_pattern = Chem.MolFromSmarts("C=CC")  # To catch double bonds flanked by carbon chains
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    num_double_bonds = len(set([(match[0], match[1]) for match in double_bond_matches]))
    if num_double_bonds != 3:
        return False, f"Found {num_double_bonds} double bonds, need exactly 3"

    # Ensure the chains are long enough, which is typical in fatty acids
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length < 12:
        return False, "Carbon chain too short for a fatty acid"

    return True, "Contains a fatty acid structure with three double bonds"