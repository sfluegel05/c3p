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

    # 2. Identify the Ketone
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone group found"
    
    for ketone_match in ketone_matches:
        carbonyl_carbon_index = ketone_match[0]

        # 3. Identify the Alpha Carbon(s)
        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_index)
        alpha_carbon_indices = [neighbor.GetIdx() for neighbor in carbonyl_carbon.GetNeighbors() if neighbor.GetSymbol() == "C"]

        # 4. Check for Hydroxy Group
        for alpha_carbon_index in alpha_carbon_indices:
           alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_index)
           hydroxy_group_found = False
           for neighbor in alpha_carbon.GetNeighbors():
               if neighbor.GetSymbol() == "O" and neighbor.GetTotalValence() == 2:
                   hydroxy_group_found = True
                   break

           if hydroxy_group_found:
               return True, "Hydroxy group found on an alpha carbon of a ketone"
    
    return False, "No hydroxy group found on any alpha carbon of a ketone"