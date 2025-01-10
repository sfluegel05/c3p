"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester

Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is an ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define butyrate ester pattern: linear four-carbon chain attached to carbonyl carbon in ester linkage
    butyrate_ester_pattern = Chem.MolFromSmarts("[CH3][CH2][CH2][C](=O)[O][#6]")
    if butyrate_ester_pattern is None:
        return False, "Invalid SMARTS pattern for butyrate ester"

    # Search for butyrate ester pattern in the molecule
    matches = mol.GetSubstructMatches(butyrate_ester_pattern)
    if matches:
        return True, "Contains butyrate ester group"
    else:
        return False, "Does not contain butyrate ester group"

__metadata__ = {   
    'chemical_class': {   
        'name': 'butyrate ester',
        'definition': 'Any carboxylic ester where the carboxylic acid component is butyric acid.'
    },
    'examples': [
        {
            'name': 'ethyl butyrate',
            'smiles': 'CCCC(=O)OCC',
            'result': True
        },
        {
            'name': 'isobutyl butyrate',
            'smiles': 'CCC(C)COC(=O)CCC',
            'result': True
        },
        {
            'name': 'butyl butyrate',
            'smiles': 'CCCCOC(=O)CCC',
            'result': True
        },
        {
            'name': 'methyl acetate',
            'smiles': 'CC(=O)OC',
            'result': False
        },
        {
            'name': 'tributyrin',
            'smiles': 'CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC',
            'result': True
        }
    ]
}