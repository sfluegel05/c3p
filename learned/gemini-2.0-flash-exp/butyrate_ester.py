"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyryl group (CH3CH2CH2C(=O)-) linked to an oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the butyrate substructure using SMARTS pattern,
    # looking for -[CH2]-[CH2]-[CH3]-(C(=O)O)- where the carbon next to the carbonyl
    # is only bonded to hydrogens or carbons.
    butyrate_pattern = Chem.MolFromSmarts("[CH2X4]-[CH2X4]-[CH3X4][CX3](=[OX1])[OX2]")

    matches = mol.GetSubstructMatches(butyrate_pattern)
    if not matches:
          return False, "No butyrate substructure found"

    # Check that the carbonyl carbon is only bonded to one other carbon atom (and 1 oxygen, and 1=O)
    # if so the pattern matches correctly
    for match in matches:
        carbonyl_c_idx = match[3]
        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        
        # Check number of C-bonds to carbonyl carbon
        carbon_bonds = 0
        for neighbor in carbonyl_c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_bonds+=1
        
        if carbon_bonds != 3:
            return False, "Carbonyl carbon does not have 3 C-bonds"
            
    # Check that no C=O group is directly attached to carbonyl carbon
    
    for match in matches:
         carbonyl_c_idx = match[3]
         carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
         for neighbor in carbonyl_c_atom.GetNeighbors():
              if neighbor.GetAtomicNum() == 6:
                   for neighbor2 in neighbor.GetNeighbors():
                         if neighbor2.GetAtomicNum() == 8 and neighbor2.GetBondBetweenAtoms(neighbor.GetIdx(), neighbor2.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                              return False, "Butyrate carbon adjacent to another C=O"
    return True, "Contains butyrate ester substructure"