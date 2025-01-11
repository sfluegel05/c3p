"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is any aliphatic monocarboxylic acid typically with an even-numbered carbon chain
    that may be saturated or unsaturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check the length of the carbon chain excluding the carboxyl group
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    
    # Assuming that one carbon is in the carboxyl group; the rest are the aliphatic chain
    aliphatic_chain_length = carbon_count - 1
    if aliphatic_chain_length < 4 or aliphatic_chain_length > 28:
        return False, f"Aliphatic chain length is {aliphatic_chain_length}, not within the 4 to 28 range"

    # Confirm the major backbone is composed of carbon and hydrogen
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (1, 6, 8):
            return False, f"Contains non-carbon and non-hydrogen atoms: atomic number {atomic_num}"

    return True, "Contains a carboxylic acid group and a suitable aliphatic hydrocarbon chain characteristic of fatty acids"