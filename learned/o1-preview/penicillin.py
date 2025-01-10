"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as any substituted penam containing two methyl groups at position 2,
    a carboxylate group at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define penam core pattern (beta-lactam fused to thiazolidine ring)
    penam_core_smarts = '[N;R1]1C(=O)[C@@H]2SC(C)(C)[C@@H](N2C1=O)'
    penam_core = Chem.MolFromSmarts(penam_core_smarts)
    if not mol.HasSubstructMatch(penam_core):
        return False, "Penam core not found"

    # Check for two methyl groups at position 2
    # Position 2 is the carbon next to the sulfur in the thiazolidine ring
    methyl_groups = Chem.MolFromSmarts('[C@@H](C)(C)[S]')
    methyl_matches = mol.GetSubstructMatches(methyl_groups)
    if not methyl_matches:
        return False, "Two methyl groups at position 2 not found"

    # Check for carboxylate group at position 3
    # Position 3 is the carbon next to position 2 in the thiazolidine ring
    carboxylate_smarts = '[C@H](N1C(=O)[C@@H]2SC(C)(C)[C@@H]2N1=O)[C](=O)[O-,O]'
    carboxylate_group = Chem.MolFromSmarts(carboxylate_smarts)
    if not mol.HasSubstructMatch(carboxylate_group):
        return False, "Carboxylate group at position 3 not found"

    # Check for carboxamido group at position 6
    # Position 6 is the carbonyl carbon in the beta-lactam ring attached to the side chain
    carboxamido_smarts = '[C@@H]1C(=O)N2[C@@H](C(O))C(C)(C)S[C@]12[C][(N)C(=O)]'
    carboxamido_group = Chem.MolFromSmarts('NC(=O)[C@@H]1N2C(=O)[C@H]2SC(C)(C)[C@H]1O')
    if not mol.HasSubstructMatch(carboxamido_group):
        return False, "Carboxamido group at position 6 not found"

    return True, "Molecule matches penicillin structure with required substituents"

__metadata__ = {
    'chemical_class': {
        'name': 'penicillin',
        'definition': 'Any member of the group of substituted penams containing two methyl substituents at position 2, a carboxylate substituent at position 3 and a carboxamido group at position 6.'
    }
}