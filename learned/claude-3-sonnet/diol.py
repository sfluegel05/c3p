"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23553 diol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound containing exactly two hydroxyl (-OH) groups.

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

    # Count total number of OH groups
    # [OX2H1] matches oxygen with 2 connections, 1 hydrogen - i.e. hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    total_oh = len(oh_matches)

    if total_oh != 2:
        return False, f"Found {total_oh} hydroxyl groups, need exactly 2"

    # Get the atoms that the OH groups are attached to
    oh_carbons = []
    for match in oh_matches:
        oh_oxygen = mol.GetAtomWithIdx(match[0])
        for neighbor in oh_oxygen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                oh_carbons.append(neighbor)

    # Check if any OH is part of a carboxylic acid
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxylic acid OH group(s)"

    # Check if any OH is part of an enol
    enol_pattern = Chem.MolFromSmarts("[CX3]=[CX3][OX2H1]")
    if mol.HasSubstructMatch(enol_pattern):
        return False, "Contains enol OH group(s)"

    # Verify both OH groups are attached to carbons
    if len(oh_carbons) != 2:
        return False, "Not all hydroxyl groups are attached to carbon atoms"

    # Success - determine if it's an aliphatic diol (glycol) or not
    is_aliphatic = all(not atom.GetIsAromatic() for atom in oh_carbons)
    if is_aliphatic:
        return True, "Aliphatic diol (glycol) with two hydroxyl groups"
    else:
        return True, "Diol with two hydroxyl groups, at least one on aromatic carbon"