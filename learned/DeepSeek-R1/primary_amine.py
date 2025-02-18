"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:15734 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine has exactly one hydrocarbyl group attached to an NH2 group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for primary amine: NH2 group attached to a carbon
    # [N;H2;D1] ensures nitrogen has two hydrogens and one bond
    # connected to a carbon (not part of amide, etc.)
    primary_amine_pattern = Chem.MolFromSmarts("[N;H2;D1;+0]-[#6]")
    matches = mol.GetSubstructMatches(primary_amine_pattern)
    
    if len(matches) > 0:
        return True, f"Found {len(matches)} primary amine group(s)"
    else:
        return False, "No primary amine group (NH2 with single carbon bond) found"