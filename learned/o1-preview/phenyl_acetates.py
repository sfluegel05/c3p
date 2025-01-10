"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester obtained by formal condensation of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phenyl acetate ester linkage
    # Aromatic carbon connected to oxygen connected to acetyl group
    phenyl_acetate_pattern = Chem.MolFromSmarts('[c][O][C](=O)C')

    # Search for phenyl acetate ester linkage in the molecule
    matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
    if matches:
        return True, "Contains phenyl acetate ester linkage"
    else:
        return False, "No phenyl acetate ester linkage found"

__metadata__ = {
    'chemical_class': {
        'name': 'phenyl acetates',
        'definition': 'An acetate ester obtained by formal condensation of the carboxy group of acetic acid with the hydroxy group of any phenol.'
    }
}