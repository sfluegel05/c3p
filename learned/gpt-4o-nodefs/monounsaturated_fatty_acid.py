"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monounsaturated fatty acid, False otherwise.
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for terminal carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for exactly one isolated carbon-carbon double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 1:
        return False, f"Expected 1 isolated carbon-carbon double bond, found {len(double_bond_matches)}"

    # Detection of ring structures
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Unexpected ring structures found"

    # Ensure carbon chain is linear and has minimal branching
    carbon_chain = Chem.MolFromSmarts("[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No linear carbon chain detected"

    # Ensure at least six carbons (as inferred from examples including shorter chains)
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count < 6:
        return False, f"Carbon chain too short for typical fatty acid, found {carbon_count} carbons"
    
    return True, "Contains the correct structure for a monounsaturated fatty acid"