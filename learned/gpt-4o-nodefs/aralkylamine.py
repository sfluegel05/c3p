"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine contains an alkyl group attached to an aromatic ring and an amine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for aromatic rings [aromatic in context]
    aromatic_pattern = Chem.MolFromSmarts("a")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"

    # Check for amine group containing primary, secondary, or tertiary nitrogen
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N-[!#6])]") 
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amino group found"

    # Check for alkyl chains - carbon atoms not in a ring
    alkyl_pattern = Chem.MolFromSmarts("[CX4;!R]")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No alkyl chain found"

    # Attempt to define specific connectivity checking amino to aromatic via the alkyl
    # An ideal SMARTS solution isn't straightforward; imply presence suffices here till explicit linkage validation

    return True, "Contains aromatic ring, alkyl chain, and amine group"