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

    # Check for a long aliphatic chain (at least 12 carbons)
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if len(c_atoms) < 12:
        return False, f"Aliphatic chain too short ({len(c_atoms)} carbons)"

    # Check for amino group attached to carbon
    amino_pattern = Chem.MolFromSmarts("[CX4][NX3;H2]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No primary amino group attached to carbon found"

    # Check for hydroxyl groups attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups attached to carbon found"

    # Check for amino alcohol sequence (C-OH-C-N)
    amino_alcohol_pattern = Chem.MolFromSmarts("[C;!R][C;!R](O)[C;!R](N)")
    if not mol.HasSubstructMatch(amino_alcohol_pattern):
        return False, "No amino alcohol sequence (C-OH-C-N) found"

    # Allow for unsaturations and additional hydroxyl groups
    # Do not enforce saturation of the aliphatic chain
    # Additional hydroxyl groups are acceptable

    # If all checks pass, classify as sphingoid
    return True, "Molecule matches sphingoid structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'sphingoid',
        'definition': 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.'
    }
}