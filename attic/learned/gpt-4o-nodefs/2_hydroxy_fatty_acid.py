"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid contains a carboxylic acid group and a hydroxyl group on the second carbon of a hydrocarbon chain.

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

    # Look for the carboxylic acid group pattern -COOH
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for the hydroxyl group on the second carbon
    hydroxyl_on_second_carbon_pattern = Chem.MolFromSmarts("C[CH](O)C(=O)[OH1]")
    if not mol.HasSubstructMatch(hydroxyl_on_second_carbon_pattern):
        return False, "No hydroxyl group on the second carbon"

    # Assess if chain is sufficiently long
    chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1

    if chain_length < 6:  # Typically, fatty acids have at least 6 carbons
        return False, "Insufficient carbon chain length for a fatty acid"

    return True, "Contains a carboxylic acid group with a hydroxyl group on the second carbon"