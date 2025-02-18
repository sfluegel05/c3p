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

    # Improved SMARTS pattern: ketone with adjacent carbon containing hydroxyl (O with exactly 1 H)
    # [CX3]=O (ketone) adjacent to a carbon connected to [O;H1] (hydroxyl)
    pattern = Chem.MolFromSmarts('[CX3](=O)-[C;!$(C=O)]([O;H1])')
    if not pattern:
        return False, "Invalid SMARTS pattern"

    # Find all matches of the pattern
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No ketone with alpha-hydroxy group detected"

    # Check if any matched ketone is valid (not part of acid/ester/amide)
    for match in matches:
        ketone_carbon = mol.GetAtomWithIdx(match[0])
        # Verify ketone environment (should have 2 carbons attached)
        carbon_neighbors = sum(1 for n in ketone_carbon.GetNeighbors() if n.GetAtomicNum() == 6)
        if carbon_neighbors >= 2:
            # Check alpha carbon's hydroxyl isn't part of conjugated system
            alpha_carbon = mol.GetAtomWithIdx(match[1])
            for bond in alpha_carbon.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(alpha_carbon).GetAtomicNum() == 8:
                    continue  # skip conjugated cases
            return True, "Contains ketone with hydroxyl group on adjacent alpha-carbon"

    return False, "No valid ketone with alpha-hydroxy group detected"