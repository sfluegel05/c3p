"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is a phenol with a methoxy group at the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for a benzene ring with an OH and OCH3
    guaiacol_base_pattern = Chem.MolFromSmarts("c1ccccc1[OH]")
    methoxy_pattern = Chem.MolFromSmarts("c1ccccc1[OCH3]")
    
    if guaiacol_base_pattern is None or methoxy_pattern is None:
       return None, "Invalid SMARTS pattern"

    # Check for the presence of both -OH and -OCH3
    if not mol.HasSubstructMatch(guaiacol_base_pattern) or not mol.HasSubstructMatch(methoxy_pattern):
         return False, "Does not have both a phenol and a methoxy group"

    # Get the indices of all atoms in the mol
    oh_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[c,C][OH]"))
    och3_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[c,C][OCH3]"))
    
    # Check for ortho relationship.
    # If the indices match, it means the OH and OCH3 are on the same ring, and therefore can be evaluated for ortho
    for oh_match in oh_matches:
        for och3_match in och3_matches:
           if mol.GetAtomWithIdx(oh_match[0]).GetOwningMol()==mol.GetAtomWithIdx(och3_match[0]).GetOwningMol():
               oh_atom_index = oh_match[0]
               och3_atom_index = och3_match[0]

               #Get neighbors of the OH attached atom
               oh_neighbors = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(oh_atom_index).GetNeighbors()]
               # Check if the methoxy carbon is one of the oh neighbours
               if och3_atom_index in oh_neighbors:
                 return True, "Has Guaiacol substructure"

               #Otherwise check the neighbours of the methoxy carbon if the oh carbon is there
               och3_neighbors = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(och3_atom_index).GetNeighbors()]
               if oh_atom_index in och3_neighbors:
                return True, "Has Guaiacol substructure"


    return False, "Does not have Guaiacol substructure"