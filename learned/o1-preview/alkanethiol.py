"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:29255 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define alkanethiol SMARTS pattern: sulfur bonded to hydrogen and sp3 carbon
    alkanethiol_pattern = Chem.MolFromSmarts("[#16H1][CX4H,CH2,CH]")
    if not mol.HasSubstructMatch(alkanethiol_pattern):
        return False, "No -SH group attached to an alkyl group found"

    # Define amide bond SMARTS pattern: N-C(=O)
    amide_pattern = Chem.MolFromSmarts("N-C(=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide bond, not an alkanethiol"

    # Optional: Exclude other functional groups
    # For example, exclude carboxylic acids, esters, nitriles, etc.
    excluded_patterns = [
        Chem.MolFromSmarts("C(=O)[OH]"),  # carboxylic acid
        Chem.MolFromSmarts("C(=O)O"),     # ester
        Chem.MolFromSmarts("C#N"),        # nitrile
        Chem.MolFromSmarts("S(=O)(=O)"),  # sulfonyl
        Chem.MolFromSmarts("N[CX3]=[CX3]"),  # imine
        Chem.MolFromSmarts("P"),          # phosphorous-containing groups
        Chem.MolFromSmarts("[!#6;!#1]"),  # any atom that is not carbon or hydrogen
    ]
    for pattern in excluded_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains excluded functional groups"

    return True, "Contains an -SH group attached to an alkyl group"