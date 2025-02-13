"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5beta-cholanic acid occurring in bile,
    typically with -OH groups and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core bile acid skeleton with 5β-configuration
    core_pattern = Chem.MolFromSmarts('[C@H]12C[C@H](O)C3(C)[CH2]C4CCC5(C)C=CCC5(C)C4(C)C3CC1O2')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No consistent 5β-cholanic acid backbone detected"

    # Presence of hydroxy groups
    hydroxy_group = Chem.MolFromSmarts('[CX4][OH]')
    ho_count = len(mol.GetSubstructMatches(hydroxy_group))
    if ho_count < 1:
        return False, "No sufficient hydroxy groups identified"

    # Presence of a carboxylic acid group or its derivatives
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OX1H0-,OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or derivative group found"

    # Analyzing possible amide linkage indicating glycine or taurine conjugation
    amide_pattern = Chem.MolFromSmarts('N-C(=O)')
    if mol.HasSubstructMatch(amide_pattern):
        return True, "Detected bile acid with an amide linkage typical of glycine or taurine conjugation"

    return True, "Contains hydroxy-5β-cholanic acid structure with required functional groups"