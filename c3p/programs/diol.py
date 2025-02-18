"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    
    # Find all matches for the hydroxyl pattern
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Count the number of hydroxyl groups
    hydroxyl_count = len(hydroxyl_matches)

    # Check if there are exactly two hydroxyl groups
    if hydroxyl_count != 2:
        return False, f"Contains {hydroxyl_count} hydroxyl groups, need exactly 2"

    # Ensure that the hydroxyl groups are not part of carboxylic acids or esters
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")

    # Check if any hydroxyl group is part of a carboxylic acid or ester
    for match in hydroxyl_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    for bond in neighbor.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            if bond.GetBeginAtom().GetAtomicNum() == 8 or bond.GetEndAtom().GetAtomicNum() == 8:  # Oxygen
                                return False, "Hydroxyl groups are part of carboxylic acids or esters"

    return True, "Contains exactly two hydroxyl groups and no disqualifying functional groups"