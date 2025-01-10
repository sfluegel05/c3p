"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose with the 'D' configuration at position 5.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly 6 carbons in the molecule for hexoses
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, "Not a hexose, incorrect number of carbon atoms"
    
    # Define common SMARTS patterns for D-hexoses in cyclic forms considering stereochemistry
    d_hexose_pyranose_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O')  # D-pyranose form
    d_hexose_furanose_pattern = Chem.MolFromSmarts('O1[C@@H]([C@H](O)[C@@H](O)C1O)[C@H](O)CO')  # D-furanose form
    d_hexose_open_chain_pattern = Chem.MolFromSmarts('[H]C(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO')  # Open-chain

    # Check for D-hexose structure matches
    if mol.HasSubstructMatch(d_hexose_pyranose_pattern) or mol.HasSubstructMatch(d_hexose_furanose_pattern) or mol.HasSubstructMatch(d_hexose_open_chain_pattern):
        return True, "Molecule identified as D-hexose with stereochemistry at carbon-5"
    
    return False, "No D-hexose structure found or incorrect stereochemistry"