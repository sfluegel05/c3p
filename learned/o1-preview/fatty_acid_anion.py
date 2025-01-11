"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation
    of the carboxylic acid group of the corresponding fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is an anion (negative charge)
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if total_charge >= 0:
        return False, "Molecule is not an anion"

    # Look for carboxylate group(s) (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "No deprotonated carboxylate group found"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Set minimum carbon count for fatty acids (e.g., 4 carbons)
    if c_count < 4:
        return False, f"Molecule has only {c_count} carbons, too few to be a fatty acid anion"

    # Calculate ratio of heteroatoms (exclude C and H)
    total_atoms = mol.GetNumAtoms()
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1,6]]
    heteroatom_ratio = len(heteroatoms) / total_atoms

    # Allow some heteroatoms for functional groups (e.g., hydroxyl groups)
    if heteroatom_ratio > 0.3:
        return False, f"Too many heteroatoms ({len(heteroatoms)} out of {total_atoms} atoms)"

    return True, "Molecule is a fatty acid anion"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58945',
                              'name': 'fatty acid anion',
                              'definition': 'The conjugate base of a fatty acid, arising '
                                            'from deprotonation of the carboxylic acid group '
                                            'of the corresponding fatty acid.',
                              'parents': ['CHEBI:25695']}}