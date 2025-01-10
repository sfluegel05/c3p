"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    Prenylquinones have a quinone core structure with prenyl chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quinone core
    quinone_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone core structure found"

    # Define prenyl side chain pattern (isoprene unit)
    prenyl_pattern = Chem.MolFromSmarts("C=C(C)C")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if len(prenyl_matches) < 1:
        return False, "No prenyl-like side chains found"

    return True, "Contains quinone core with prenyl side chain(s)"

# Example usage
example_smiles = "COC1=CC(=O)C=C(C\C=C(/C)CC\C=C(/C)CCC=C(C)C)C1=O"
is_prenylquinone(example_smiles)