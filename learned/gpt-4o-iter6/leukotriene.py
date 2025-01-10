"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes have a characteristic C20 polyunsaturated backbone with at
    least four double bonds, predominantly conjugated, and usually terminate
    with a carboxylic acid group.

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

    # Count carbon atoms - allow for derivatives with extra groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Only {c_count} carbon atoms, expected at least 20 including substitutions for leukotrienes"

    # Identify double bonds, allowing flexibility in derivatives
    db_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if db_count < 4:
        return False, f"Contains {db_count} double bonds; expected at least 4 for polyunsaturation characteristic of leukotrienes"

    # Check for three or more conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts('C=,C=,C=C=C')
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Four non-complex-conjugated double bonds required; no valid pattern found"

    # Check presence of a terminal carboxylic acid chain considering typical length
    terminal_acid_pattern = Chem.MolFromSmarts('CCCC(O)=O')
    if not mol.HasSubstructMatch(terminal_acid_pattern):
        return False, "Lack of typical terminal carboxylic acid as found in leukotriene structures"

    return True, "Matches leukotriene profile: C20 unsaturated backbone with conjugated double bonds and terminal functionality"