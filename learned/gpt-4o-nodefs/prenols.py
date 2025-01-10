"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols generally consist of multiple isoprene units and contain an alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Missing alcohol group"

    # Check for isoprene units. Isoprene units in prenols are often arranged with repeated C(=C)-C-C=C motifs
    isoprene_pattern = Chem.MolFromSmarts('C(=C)-C-C=C')
    isoprene_count = len(mol.GetSubstructMatches(isoprene_pattern))
    if isoprene_count < 1:
        return False, f"Insufficient isoprene-like units (found {isoprene_count})"

    return True, "Contains alcohol group and isoprene units, typical of prenols"