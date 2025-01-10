"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Only {c_count} carbon atoms, minimum expected is 20 for leukotrienes"

    # Identify double bonds
    db_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if db_count < 4:
        return False, f"Contains {db_count} double bonds; expected at least 4 for polyunsaturation characteristic"

    # Check for three conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    quoted = Chem.MolFromSmarts("C=C-C=C-C=C") if conjugated_pattern else None
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 1:
        return False, "Requires three or more conjugated double bonds; none found"

    # Check presence of a terminal carboxylic acid or similar functionality
    # Expanded pattern for variations of terminal group
    terminal_acid_pattern = Chem.MolFromSmarts("[C,c](=O)[O,o]")
    if not mol.HasSubstructMatch(terminal_acid_pattern):
        return False, "Lack of typical terminal carboxylic acid or variant functionality found in leukotriene structures"

    return True, "Matches leukotriene profile: C20 unsaturated backbone with conjugated double bonds and terminal functionality"