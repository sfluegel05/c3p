"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: alpha-hydroxy ketone (CHEBI:59805)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone has a ketone group with a hydroxyl group on the adjacent (alpha) carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all ketone groups (any C=O)
    ketone_pattern = MolFromSmarts('[CX3]=O')
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone group found"

    # Check each ketone's adjacent carbons for hydroxyl groups
    for match in ketone_matches:
        ketone_carbon = mol.GetAtomWithIdx(match[0])
        # Get neighboring carbons (alpha carbons)
        alpha_carbons = [neighbor for neighbor in ketone_carbon.GetNeighbors() 
                        if neighbor.GetAtomicNum() == 6]
        
        for alpha in alpha_carbons:
            # Check for hydroxyl group (-OH) on alpha carbon
            hydroxyl_found = False
            for bond in alpha.GetBonds():
                other = bond.GetOtherAtom(alpha)
                # Verify it's an -OH group (O with single bond and at least one H)
                if (other.GetAtomicNum() == 8 and 
                    bond.GetBondType() == Chem.BondType.SINGLE and 
                    other.GetTotalNumHs() >= 1):
                    hydroxyl_found = True
                    break
            if hydroxyl_found:
                return True, (f"Hydroxyl group found on alpha carbon adjacent to "
                             f"ketone at position {ketone_carbon.GetIdx()+1}")

    return False, "No hydroxyl group on alpha carbon of any ketone"