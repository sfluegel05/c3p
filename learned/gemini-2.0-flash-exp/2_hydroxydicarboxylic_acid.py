"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid with a hydroxyl group on the alpha carbon
    relative to *each* carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True and a reason string if the molecule is a 2-hydroxydicarboxylic acid,
                         False and a reason string otherwise, or None, None if input is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # 1. Check for two carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 2:
       return False, f"Molecule does not have 2 carboxyl groups. Found {len(carboxyl_matches)}"


    # 2. Check for at least one hydroxyl on an alpha carbon relative to *each* carboxyl
    alpha_hydroxyl_count = 0
    for carboxyl_match in carboxyl_matches:
        carboxyl_carbon_idx = carboxyl_match[0] # Get index of the carboxyl carbon
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        alpha_carbon_idx = -1
        for neighbor in carboxyl_carbon.GetNeighbors():
              #Find the alpha carbon, which should be a carbon and not the oxygen of the carboxyl
              if neighbor.GetAtomicNum() == 6:
                 alpha_carbon_idx = neighbor.GetIdx()
                 break
        if alpha_carbon_idx == -1:
           return False, "Carboxyl group does not have an alpha carbon"

        # Check if this specific alpha carbon has a hydroxyl group
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

        # Using a smarter SMARTS that checks for alcohol carbons -[C](-O)-
        alcohol_pattern = Chem.MolFromSmarts("[CX4]-[OX2H1]") #  Explicitly look for C-O-H 
        
        for hydroxyl_match in mol.GetSubstructMatches(alcohol_pattern):
              hydroxyl_carbon_idx = mol.GetAtomWithIdx(hydroxyl_match[0]).GetIdx() # index of the alcohol carbon

              if hydroxyl_carbon_idx == alpha_carbon_idx: # if the carbon attached to the alcohol IS the alpha carbon of the carboxyl
                alpha_hydroxyl_count+=1
                break
            
    if alpha_hydroxyl_count != 2:
      return False, f"Not all alpha carbons have a hydroxyl. Found {alpha_hydroxyl_count}"


    #4 Success
    return True, "Molecule is a 2-hydroxydicarboxylic acid"