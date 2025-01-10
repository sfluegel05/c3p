"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone typically has a 2,3-dimethoxybenzoquinone core with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the 2,3-dimethoxy-5-methylbenzoquinone core
    benzoquinone_pattern = Chem.MolFromSmarts("COc1c(C)cc(=O)c(=O)c1OC")  # Ensure methoxy groups and methyl position
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No 2,3-dimethoxybenzoquinone core found"
    
    # Look for polyprenoid side chain at position 6 (identify by repeating isoprene units, allowing isomerism)
    isoprene_unit = Chem.MolFromSmarts("CC(C)=C")  # Typical isoprene pattern
    matches = mol.GetSubstructMatches(isoprene_unit)
    if not matches or len(matches) < 1:  # At least one isoprene unit expected in side chain
        return False, "Insufficient polyprenoid side chain found"
    
    # Look for quinoid functional group to support ubiquinone stability
    quinoid_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")  # Basic quinoid pattern check
    if not mol.HasSubstructMatch(quinoid_pattern):
        return False, "No redox-active quinoid group detected"

    return True, "Contains ubiquinone core with reasonable polyprenoid side chain and quinoid group"