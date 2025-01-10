"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is classified as a wax based on its SMILES string.
    Waxes are characterized by long-chain fatty acids esterified to long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester bond pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    # Check for presence of ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond present"

    # Count total number of carbon atoms
    num_of_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check if the total carbon atom count is relatively high, suggesting a wax
    if num_of_carbons >= 20:
        return True, "Molecule with ester bond and high carbon count consistent with wax structure"

    return False, "Ester bond found but total carbon count too low to be typical of a wax"