"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule belongs to the acrovestone class based on its SMILES string.
    Acrovestone compounds are polyphenols isolated from Acronychia pedunculata, typically characterized by an isoflavone core with glycosylated substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to acrovestone class, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general isoflavone core SMARTS pattern (3-phenylchromen-4-one)
    isoflavone_pattern = Chem.MolFromSmarts('c1ccc(cc1)c2cc(=O)c3ccccc3o2')
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Check for glycosylation (sugar moieties attached via oxygen)
    # Look for O-glycosidic bond: aromatic carbon - oxygen - sugar ring
    glycoside_pattern = Chem.MolFromSmarts('c[O][C;R]')
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosylation detected"

    # Additional checks can include counting hydroxyl and methoxy groups if needed
    # For now, assume detection of isoflavone core and glycosylation is sufficient

    reason = "Molecule matches acrovestone class (isoflavone core with glycosylation)"
    details = []
    details.append("Isoflavone core present")
    details.append("Glycosylation detected")

    return True, reason + "; " + "; ".join(details)

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'acrovestone',
        'definition': 'A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata can be added as required
}