"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:2468 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone has a ketone group (C=O) with a hydroxyl (-OH) on the adjacent (alpha) carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised SMARTS pattern to accurately capture alpha-hydroxy ketone structure:
    # [CX3]=O (ketone) adjacent to a carbon that has a hydroxyl group directly attached
    pattern = Chem.MolFromSmarts('[CX3](=O)-[C]([OH])')
    if not pattern:
        return False, "Invalid SMARTS pattern"

    # Find all potential alpha-hydroxy ketone matches
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No ketone with alpha-hydroxy group detected"

    # Verify ketone validity (must have >=2 carbon neighbors) and check OH group type
    for match in matches:
        ketone_carbon = mol.GetAtomWithIdx(match[0])
        carbon_neighbors = sum(1 for bond in ketone_carbon.GetBonds() 
                              if bond.GetOtherAtom(ketone_carbon).GetAtomicNum() == 6)
        
        if carbon_neighbors >= 2:  # Valid ketone (not aldehyde/acid/ester)
            # Final validation: Ensure OH is not part of conjugated system/carboxylic acid
            alpha_carbon = mol.GetAtomWithIdx(match[1])
            for bond in alpha_carbon.GetBonds():
                other = bond.GetOtherAtom(alpha_carbon)
                if other.GetAtomicNum() == 8 and other.GetTotalNumHs() == 1 and bond.GetBondType() == Chem.BondType.SINGLE:
                    return True, "Contains ketone with hydroxyl group on adjacent alpha-carbon"

    return False, "No valid ketone with alpha-hydroxy group detected"