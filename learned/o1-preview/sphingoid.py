"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sphingoid backbone pattern:
    # Long aliphatic chain with amino group at C2, hydroxyls at C1 and C3
    sphingoid_pattern = Chem.MolFromSmarts("""
    [#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-
    [C@H](O)[C@H](N)[C@H](O)CO
    """)
    if sphingoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if molecule matches the sphingoid pattern
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "Molecule does not match sphingoid backbone pattern"

    # Alternatively, check for components individually
    # 1. Long aliphatic chain (at least 12 carbons)
    aliphatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 2]
    if len(aliphatic_carbons) < 12:
        return False, f"Aliphatic chain too short ({len(aliphatic_carbons)} carbons)"

    # 2. Amino group at C2 position
    amino_pattern = Chem.MolFromSmarts("[C;!R][C;!R](N)")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group at C2 position found"

    # 3. Hydroxyl groups at C1 and C3 positions
    hydroxyl_pattern = Chem.MolFromSmarts("[C;!R](O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Less than two hydroxyl groups found"

    # 4. Check for unsaturations (double bonds) and additional hydroxy groups
    # Unsaturation is allowed, so no need to enforce saturation
    # Additional hydroxyl groups are allowed

    # If all checks pass, classify as sphingoid
    return True, "Molecule matches sphingoid structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'sphingoid',
        'definition': 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.'
    }
}