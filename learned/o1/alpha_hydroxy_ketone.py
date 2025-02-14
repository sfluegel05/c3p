"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-Hydroxy Ketone
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone containing a hydroxy group on 
    the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ketone groups (exclude aldehydes, carboxylic acids, esters)
    ketone_pattern = Chem.MolFromSmarts('C(=O)C')  # Ketone carbonyl group
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    if not ketone_matches:
        return False, "No ketone group found"

    for match in ketone_matches:
        ketone_carbon_idx = match[0]  # Carbonyl carbon
        alpha_carbons = []

        # Get adjacent carbons (alpha-carbons)
        for nbr in mol.GetAtomWithIdx(ketone_carbon_idx).GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbons.append(nbr.GetIdx())

        # Check if any alpha-carbon has a hydroxy group attached
        for alpha_carb_idx in alpha_carbons:
            has_hydroxy = False
            alpha_carbon = mol.GetAtomWithIdx(alpha_carb_idx)
            for nbr in alpha_carbon.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:  # Oxygen atom with single bond (hydroxy)
                    has_hydroxy = True
                    break
            if has_hydroxy:
                return True, "Contains an alpha-hydroxy ketone group"

    return False, "Does not contain an alpha-hydroxy ketone group"