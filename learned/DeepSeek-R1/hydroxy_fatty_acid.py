"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: CHEBI hydroxy fatty acid
"""
from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid (carboxylic acid with aliphatic chain) 
    carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_matches:
        return False, "No carboxylic acid group"

    # Collect oxygen atoms in carboxylic acid groups
    carboxylic_oxygens = {match[2] for match in carboxylic_matches}  # Index of -OH oxygen in COOH

    # Check for hydroxyl groups not part of carboxylic acid
    hydroxyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:  # Hydroxyl oxygen
            if atom.GetIdx() not in carboxylic_oxygens:
                hydroxyl_found = True
                break

    if not hydroxyl_found:
        return False, "No hydroxyl groups outside carboxylic acid"

    return True, "Carboxylic acid with hydroxyl substituent(s)"