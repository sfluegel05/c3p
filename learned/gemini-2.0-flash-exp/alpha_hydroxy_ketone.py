"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone with a hydroxy group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Define SMARTS pattern for ketone and neighboring carbon
    ketone_alpha_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")
    
    # 3. Get all matches of the ketone and neighbor.
    matches = mol.GetSubstructMatches(ketone_alpha_pattern)

    # 4. If there is no ketone, return False.
    if not matches:
      return False, "No ketone found."
    
    # 5. Check if any of the alpha carbons have a hydroxy group.
    for match in matches:
        ketone_carbon_idx = match[0] # index of the ketone carbon
        alpha_carbon_idx = match[1] # index of the alpha carbon.
        
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

        # Get neighbors of the alpha carbon
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8: # check for oxygen
                # Check if the oxygen is an -OH (i.e. bonded to a hydrogen)
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetAtomicNum() == 1:
                        return True, "Hydroxy group found on an alpha carbon of a ketone"
    
    return False, "No hydroxy group found on any alpha carbon of a ketone"