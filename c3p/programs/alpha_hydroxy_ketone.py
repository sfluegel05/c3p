"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:2468 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts

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

    # Enhanced SMARTS pattern to capture both aliphatic and aromatic ketones
    # 1. Find ketone (C=O) with adjacent carbon (alpha)
    # 2. Alpha carbon must have hydroxyl (O with 1 H)
    # 3. Exclude conjugated systems and ester/acid environments
    pattern = MolFromSmarts('[#6]=[O]-[#6;!$(C=O)]-[O;H1]')
    if not pattern:
        return False, "Invalid SMARTS pattern"

    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No ketone with alpha-hydroxy group detected"

    # Additional validation checks for each match
    for match in matches:
        ketone_idx, alpha_idx, _, oh_idx = match
        
        # Check ketone environment - must have at least 2 carbon neighbors
        ketone_atom = mol.GetAtomWithIdx(ketone_idx)
        carbon_neighbors = sum(1 for n in ketone_atom.GetNeighbors() if n.GetAtomicNum() == 6)
        if carbon_neighbors < 2:
            continue  # Reject carbonyls in acids/esters
        
        # Check alpha carbon for conjugation (no double bonds)
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                break
        else:  # Only execute if no break (no double bonds)
            return True, "Contains ketone with hydroxyl group on adjacent alpha-carbon"

    return False, "No valid alpha-hydroxy ketone pattern found"