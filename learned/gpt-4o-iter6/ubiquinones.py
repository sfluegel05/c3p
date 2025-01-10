"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are benzoquinones with methoxy groups and a polyisoprenoid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the benzoquinone core pattern
    benzoquinone_core_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")
    if not mol.HasSubstructMatch(benzoquinone_core_pattern):
        return False, "No benzoquinone moiety found"

    # Check for methoxy groups (typically 2,3-dimethoxybenzoquinone)
    methoxy_groups_pattern = Chem.MolFromSmarts("COC1=C(OC)C=CC(=O)C=C1")
    if not mol.HasSubstructMatch(methoxy_groups_pattern):
        return False, "Required methoxy groups not found on the benzoquinone moiety"
    
    # Check for long polyisoprenoid chain (more than one isoprene unit)
    # Adjusted for simplicity but captures extendable repeats
    polyisoprenoid_pattern = Chem.MolFromSmarts("C=C(C)CCC=C")
    matches = mol.GetSubstructMatches(polyisoprenoid_pattern)
    if len(matches) < 2:
        return False, "No extended polyisoprenoid chain found"

    return True, "Matches the characteristics of a ubiquinone: contains appropriate benzoquinone moiety with methoxy groups and an attached polyisoprenoid chain"