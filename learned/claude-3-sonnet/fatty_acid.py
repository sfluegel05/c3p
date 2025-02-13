"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:36978 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_valid_smiles(smiles):
    """
    Validates a SMILES string.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if the SMILES string is valid, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except ValueError:
        return False

def is_fatty_acid(smiles):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid with a chain length
    typically ranging from 4 to 28 carbons, which may be saturated or unsaturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    try:
        # Validate SMILES string
        if not is_valid_smiles(smiles):
            return False, "Invalid SMILES string"

        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)

        # Look for carboxylic acid group
        acid_pattern = Chem.MolFromSmarts("C(=O)O")
        if not mol.HasSubstructMatch(acid_pattern):
            return False, "No carboxylic acid group found"

        # Look for aliphatic chain
        chain_pattern = Chem.MolFromSmarts("C~C~C~C")
        chain_matches = mol.GetSubstructMatches(chain_pattern)
        if not chain_matches:
            return False, "No aliphatic chain found"

        # Count number of carbon atoms in the longest chain
        longest_chain_length = max(len(set([mol.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx()
                                            for match in chain_matches
                                            for idx in match])) + 1, default=0)

        # Check chain length (typically 4 to 28 carbons)
        if not (4 <= longest_chain_length <= 28):
            return False, f"Chain length ({longest_chain_length}) is outside the typical range for fatty acids (4-28 carbons)"

        # Check for cyclic structures
        if any(ring.IsAromaticRing() for ring in mol.GetRingInfo().AtomRings()):
            return False, "Cyclic structures detected, fatty acids should be acyclic"

        # Check for branching (allow for limited branching)
        branching_pattern = Chem.MolFromSmarts("CC(C)(C)C")
        branching_matches = mol.GetSubstructMatches(branching_pattern)
        if len(branching_matches) > 1:
            return False, "Excessive branching detected, fatty acids should have limited branching"

        # Additional criteria can be added here, such as checking for specific functional groups
        # or handling edge cases

        return True, "Molecule meets the structural criteria for a fatty acid"

    except Exception as e:
        return False, f"An error occurred: {str(e)}"