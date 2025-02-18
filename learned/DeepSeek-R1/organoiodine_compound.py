"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies organoiodine compounds (CHEBI:37139) - compounds containing at least one carbon-iodine bond.
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on the presence of at least one carbon-iodine bond,
    excluding cases where iodine is part of specific functional groups (e.g., sulfonic acid derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains valid C-I bond, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for any carbon-iodine bonds
    c_i_bond = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Iodine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    c_i_bond = True
                    break
            if c_i_bond:
                break

    if not c_i_bond:
        return False, "No carbon-iodine bonds found"

    # Exclude specific functional groups
    # 1. Sulfonic acid derivatives (e.g., CI attached to sulfonic acid group)
    sulfonic_pattern = Chem.MolFromSmarts("[C][I](S(=O)(=O)O)")
    if mol.HasSubstructMatch(sulfonic_pattern):
        return False, "C-I bond in sulfonic acid group"

    # 2. Check for iodo substituents on inorganic-like structures (e.g., CI directly attached to multiple oxygens)
    inorganic_pattern = Chem.MolFromSmarts("[C]([I])(O)(O)")
    if mol.HasSubstructMatch(inorganic_pattern):
        return False, "C-I bond in inorganic-like structure"

    return True, "Contains valid carbon-iodine bond"