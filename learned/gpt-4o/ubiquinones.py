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
    
    # Look for the 2,3-dimethoxybenzoquinone core (allows minor variations)
    benzoquinone_pattern = Chem.MolFromSmarts("COc1c(C)cc(=O)[#6H1]=,:[#6H1]c1=O")  # Flexible SMARTS
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No 2,3-dimethoxybenzoquinone core found"

    # Look for the polyprenoid side chain (identify by repeating isoprene units)
    isoprene_unit = Chem.MolFromSmarts("CC(C)=C")  # Typical isoprene pattern
    matches = mol.GetSubstructMatches(isoprene_unit)
    if not matches or len(matches) < 2:  # Typically expect multiple repeats
        return False, "No sufficient polyprenoid side chain found"

    # Check for redox-active quinoid group
    quinoid_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")  # Checking for quinoid pattern
    if not mol.HasSubstructMatch(quinoid_pattern):
        return False, "No redox-active quinoid group detected"

    return True, "Contains ubiquinone core with adequate polyprenoid side chain and quinoid group"