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

    # Revised pattern: Phenol with ortho methoxy group
    # Using `c1ccccc1` for benzene, ensure O and OC are ortho to each other
    guaiacol_pattern_ortho = Chem.MolFromSmarts('c1c(O)cc(OC)cc1')

    # Match for the revised ortho pattern
    if mol.HasSubstructMatch(guaiacol_pattern_ortho):
        return (True, "Molecule is classified as a guaiacol: contains phenol with ortho-position methoxy group")
    
    return (False, "Molecule lacks the specific ortho-methoxyphenol structure")