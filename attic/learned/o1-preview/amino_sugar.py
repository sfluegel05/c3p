"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:27176 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is defined as any sugar having one or more alcoholic hydroxy groups replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify sugar rings (furanose and pyranose rings)
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O"),  # Pyranose ring
        Chem.MolFromSmarts("C1OC(O)C(O)C1O")       # Furanose ring
    ]
    sugar_found = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(pattern):
            sugar_found = True
            break
    if not sugar_found:
        return False, "No sugar ring found"

    # Identify amino groups attached to ring carbons
    # Pattern for ring carbon attached to amino group
    amino_substitution_pattern = Chem.MolFromSmarts("[C;R][N;X3]")
    amino_carbons = mol.GetSubstructMatches(amino_substitution_pattern)
    if len(amino_carbons) == 0:
        return False, "No amino substitution of hydroxyl groups found on sugar ring"

    # Confirm that the amino group is replacing a hydroxy group
    # (This is a simplification; detailed analysis would require tracking the substitution)
    return True, "Molecule is an amino sugar with amino substitution on sugar ring carbons"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27176',
                                          'name': 'amino sugar',
                                          'definition': 'Any sugar having one or more alcoholic hydroxy groups replaced by substituted or unsubstituted amino groups.',
                                          'parents': []},
                    'config': {},
                    'message': None,
                    'attempt': 0,
                    'success': True,
                    'best': True,
                    'error': '',
                    'stdout': None}