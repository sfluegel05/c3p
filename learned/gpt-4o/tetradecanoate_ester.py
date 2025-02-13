"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the esterification of tetradecanoic acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define tetradecanoate ester functional group pattern
    # 13 carbon main chain with ester group at the end (-C(=O)O-)
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O")
    matches = mol.GetSubstructMatches(tetradecanoate_pattern)

    if matches:
        for match in matches:
            ester_bond = False
            # Check the connection beyond the match start for ester linkage
            ester_carbons = [mol.GetAtomWithIdx(idx) for idx in match]
            for atom in ester_carbons:
                if atom.GetAtomicNum() == 8:  # Oxygen in the ester
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in match:
                            ester_bond = True
                            break
            if ester_bond:
                return True, "Contains tetradecanoate ester group"
        return False, "Found potential tetradecanoate group, but no ester bond"
    else:
        return False, "No tetradecanoate ester group found"