"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid with a chain length of less than C6.
    If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check that the molecule is aliphatic (no aromatic rings)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if len(aromatic_atoms) > 0:
        return False, "Molecule contains aromatic rings"

    # Check for non-hydrocarbon substituents (e.g., halogens, nitrogen, sulfur, etc.)
    non_hydrocarbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8]]
    if len(non_hydrocarbon_atoms) > 0:
        return False, "Molecule contains non-hydrocarbon substituents"

    # Calculate the number of carbon atoms in the longest chain excluding the carboxylic acid carbon
    longest_chain = rdMolDescriptors.CalcLongestChain(mol)
    if longest_chain >= 6:
        return False, f"Chain length is {longest_chain}, which is too long for a short-chain fatty acid"

    # Check that the molecule is a monocarboxylic acid (only one carboxylic acid group)
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    return True, "Aliphatic monocarboxylic acid with a chain length of less than C6 and no non-hydrocarbon substituents"