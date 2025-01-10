"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: Prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quinone pattern (1,4-benzoquinone or naphthoquinone)
    quinone_pattern1 = Chem.MolFromSmarts("C1=CC(=O)C=CC1=O")
    quinone_pattern2 = Chem.MolFromSmarts("C1=CC=C(O)C(=O)C=C1")
    if not mol.HasSubstructMatch(quinone_pattern1) and not mol.HasSubstructMatch(quinone_pattern2):
        return False, "No quinone backbone found"

    # Look for isoprenoid-like motifs in the side chain (C=C-C repetition)
    prenyl_pattern = Chem.MolFromSmarts("C(=C)C")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if len(prenyl_matches) < 1:
        return False, "No prenyl side-chain detected"
    
    # (Optional) Check for additional stereochemistry or patterns specific to subclasses
    # This step can be expanded based on specific knowledge of menaquinones, phylloquinones, etc.

    return True, "Contains quinone backbone with prenyl side-chain"

# Note: This is a basic implementation to demonstrate the concept. Actual classification may require
#       more comprehensive patterns and checks to ensure accurate identification of prenylquinones.