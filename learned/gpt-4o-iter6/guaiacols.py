"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is defined as any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Look for a more specific pattern of guaiacols: a phenol with an ortho-methoxy group
    guaiacol_pattern_ortho1 = Chem.MolFromSmarts('c1cc(OC)cc(O)c1')  # Ortho methoxy to hydroxy
    guaiacol_pattern_ortho2 = Chem.MolFromSmarts('c1cc(O)cc(OC)c1')  # Ortho hydroxy to methoxy
    
    # Match for either of the ortho patterns
    if mol.HasSubstructMatch(guaiacol_pattern_ortho1) or mol.HasSubstructMatch(guaiacol_pattern_ortho2):
        return (True, "Molecule is classified as a guaiacol: contains phenol with ortho-position methoxy group")

    return (False, "Molecule lacks the specific ortho-methoxyphenol structure")