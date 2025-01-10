"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:17855 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy functional group in the alpha- or 2-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    hydroxy_alpha_generic_pattern = Chem.MolFromSmarts("CC(O)C(=O)O")  # generic pattern to match more flexible structures

    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for 2-hydroxy group in alpha position generically
    if not mol.HasSubstructMatch(hydroxy_alpha_generic_pattern):
        return False, "No 2-hydroxy group found in the alpha position"
    
    # Count carbons to ensure appropriate assessment of fatty acid chain, leniently adjusting threshold
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:  # Optionally set a lower bound threshold based on the false-negative analysis
        return False, f"Insufficient carbon chain length for a fatty acid, found {carbon_count} carbons"

    return True, "Confirmed: Contains 2-hydroxy group in the alpha position of a fatty acid chain"