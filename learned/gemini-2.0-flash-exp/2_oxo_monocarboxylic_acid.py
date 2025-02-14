"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid contains one carboxylic acid group and a carbonyl group at
    the alpha position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for monocarboxylic acid group (C(=O)[O,OH1])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,OH1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Molecule must have exactly one carboxylic acid group, found {len(carboxylic_acid_matches)}"
    
    # Get the carboxyl carbon
    carboxyl_carbon_idx = carboxylic_acid_matches[0][0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

    # Identify the alpha carbon, must be a carbon that is single bonded to the carboxylic carbon
    alpha_carbon = None
    for bond in carboxyl_carbon.GetBonds():
      if bond.GetBondType() == Chem.BondType.SINGLE:
          neighbor = bond.GetOtherAtom(carboxyl_carbon)
          if neighbor.GetSymbol() == 'C':
                alpha_carbon = neighbor
                break
    if alpha_carbon is None:
      return False, "Carboxylic acid group is not connected to another carbon"
    

    # Check if alpha carbon has a carbonyl group.
    has_carbonyl = False
    for neighbor in alpha_carbon.GetNeighbors():
        if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
            has_carbonyl = True
            break
    
    if not has_carbonyl:
       return False, "No 2-oxo group found on the alpha carbon"


    return True, "Molecule has a 2-oxo substituent and one carboxylic acid group"