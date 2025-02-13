"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are typically characterized by a 5β-steroidal backbone with hydroxyl and carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and a reason if the molecule is a bile acid, False otherwise with a reason
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broad steroidal backbone pattern for 5β-cholanic acid
    steroid_pattern = Chem.MolFromSmarts('[C@H]1(CC[C@@H]2[C@@H]3(C)CC[C@H]4[C@]3(CC[C@@H](C4)C2)C1)[C@@H](C)CCC(=O)O')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No consistent steroidal 5β-backbone detected"

    # Check for hydroxy groups
    hydroxy_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[CX4][OH]'))
    if len(hydroxy_matches) < 1:
        return False, "Not enough hydroxy groups identified"

    # Check for a carboxylic acid group or its esters/amides
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OX1H0-,OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or derivative group found"

    return True, "Contains features typical of a bile acid with 5β-steroidal structure and functional groups"