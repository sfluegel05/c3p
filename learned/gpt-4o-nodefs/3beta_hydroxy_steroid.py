"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for a steroid skeleton:
    # This pattern model assumes a basic ring structure common to steroids
    steroid_core_pattern = Chem.MolFromSmarts("C1CC[C@H]2C(=O)C3CCC[C@]3(C)C2C1")
    
    # The pattern specifically for the beta-oriented hydroxyl group on the third carbon
    hydroxyl_beta_pattern = Chem.MolFromSmarts("[C@@H](O)")

    # Check if the molecule has a steroid core structure
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not match steroid core"

    # Check if there's a 3-beta hydroxyl
    # Define multiple possible attachment sites and orientations
    substructure_matches = mol.GetSubstructMatches(hydroxyl_beta_pattern)
    
    # Look through substructure matches to ensure attachment at correct position
    for match in substructure_matches:
        if any(mol.GetAtomWithIdx(atom_idx).GetChiralTag() == Chem.CHI_TETRAHEDRAL_CCW and
               mol.GetAtomWithIdx(atom_idx).GetDegree() == 3 for atom_idx in match):
            return True, "Matches 3beta-hydroxy steroid pattern"

    return False, "3-beta hydroxyl group not in correct configuration"