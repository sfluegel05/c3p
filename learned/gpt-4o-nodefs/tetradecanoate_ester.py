"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester contains the tetradecanoate moiety, which is a 14-carbon saturated fatty acid ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a tetradecanoate ester (14C chain ending with ester linkage)
    # This pattern allows for ester linkages with various functional moieties
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)[O;!H0]")
    
    if mol.HasSubstructMatch(tetradecanoate_pattern):
        return True, "Contains tetradecanoate moiety"

    # Additional properties checks (e.g. molecular weight considerations or number of carbon atoms)
    # can be added here if needed for differentiation or validation.

    return False, "Does not contain tetradecanoate moiety"