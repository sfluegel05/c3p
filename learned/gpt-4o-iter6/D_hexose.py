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
    
    # SMARTS patterns for D-pyranoses with stereochemistry at position 5
    d_pyranose_patterns = [
        Chem.MolFromSmarts('C1O[C@H]([C@@H](O)[C@H](O)[C@H](O[C@H]1)O)'),  # Simplified generic D-pyranose pattern
    ]
    
    # SMARTS patterns for D-furanoses with stereochemistry at position 5
    d_furanose_patterns = [
        Chem.MolFromSmarts('C1O[C@@H]([C@@H](O)[C@H](O)[C@@H](O1)CO)'),  # Generic D-furanose pattern
    ]
    
    # SMARTS for open-chain D-hexoses with stereochemistry at position 5
    open_chain_pattern = Chem.MolFromSmarts('C(=O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO')

    # Check for matches to any D-hexose conformation
    for pattern in d_pyranose_patterns + d_furanose_patterns + [open_chain_pattern]:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule identified as D-hexose with stereochemistry at carbon-5"

    # Additional check for phosphates or other derivative groups that may confuse classification
    if mol.HasSubstructMatch(Chem.MolFromSmarts('OP(O)(O)=O')):  # Phosphate group
        return False, "Molecule has phosphate group, likely a derivative not a standard hexose"

    return False, "No valid D-hexose structure found, or incorrect stereochemistry"