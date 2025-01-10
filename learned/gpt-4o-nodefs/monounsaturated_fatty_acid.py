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

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for exactly one carbon-carbon double bond
    double_bond_pattern = Chem.MolFromSmarts("C=CC(=O)O")
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches < 1:
        return False, f"Failed to find carbon-carbon double bond adjacent to carboxylic acid"

    # Ensure appropriate length of carbon chain
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count < 5:
        return False, f"Carbon chain too short, found {carbon_count} carbons"
    
    # Ensure system does not misinterpret configurational idiosyncrasies
    is_cis = 'cis' in smiles or '/' in smiles
    is_trans = 'trans' in smiles or '\\' in smiles
    if is_cis and is_trans:
        return False, "SMILES implies inconsistent geometry in double bond"

    return True, "Contains the correct structure for a monounsaturated fatty acid"