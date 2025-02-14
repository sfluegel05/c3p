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
    Acrovestone compounds are polyphenols isolated from Acronychia pedunculata, with an isoflavone core and glycosylated substituents.

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

    # Define a more accurate isoflavone core SMARTS pattern (3-phenylchromen-4-one)
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1c3ccccc3')
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Check for glycosylation (sugar moieties)
    # Using a pattern for glycosidic linkage (oxygen linking sugar to core)
    glycoside_pattern = Chem.MolFromSmarts('[#6]-[#8]-[C@@H]1[C@@H]([#8H])[C@H]([#8H])[C@@H]([#8H])[C@H]1[#8H]')
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosylation detected"

    # Count hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Count methoxy groups (-OCH3)
    methoxy_pattern = Chem.MolFromSmarts('[OX2][CH3]')
    num_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))

    # Count phenolic hydroxyl groups (hydroxyl attached to aromatic ring)
    phenol_pattern = Chem.MolFromSmarts('c[OX2H]')
    num_phenol = len(mol.GetSubstructMatches(phenol_pattern))

    reason = "Molecule matches acrovestone class (isoflavone core with glycosylation)"
    details = []
    details.append("Isoflavone core present")
    details.append("Glycosylation detected")
    details.append(f"Number of hydroxyl groups: {num_hydroxyl}")
    details.append(f"Number of methoxy groups: {num_methoxy}")
    details.append(f"Number of phenolic hydroxyl groups: {num_phenol}")

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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata can be added as required
}