"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    Azoles are monocyclic 5-membered heterocycles containing at least one nitrogen and
    possibly other heteroatoms (N, S, O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Sanitize SMILES string
    try:
      mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except:
      return False, "Invalid SMILES string"

    if mol is None:
        return False, "Invalid SMILES string"

    # Check if it's monocyclic
    if Chem.GetSSSR(mol) != 1:
        return False, "Molecule is not monocyclic"
    
    # Find 5-membered rings, each atom in the ring
    azole_pattern = Chem.MolFromSmarts("[R5][R5][R5][R5][R5]")
    if not mol.HasSubstructMatch(azole_pattern):
          return False, "No 5-membered ring found"

    #check for at least one nitrogen atom
    nitrogen_pattern = Chem.MolFromSmarts("[n;R5]")
    if not mol.HasSubstructMatch(nitrogen_pattern):
        return False, "No nitrogen atom in the 5-membered ring"

    #Check for no more than 2 other heteroatoms (N, O, S). Max 3 heteroatoms, min 1 nitrogen.
    other_hetero_pattern = Chem.MolFromSmarts("[o,s;R5]")
    other_hetero_matches = mol.GetSubstructMatches(other_hetero_pattern)
    if len(other_hetero_matches) > 2:
        return False, "Too many heteroatoms in the ring"
        

    return True, "Contains a monocyclic 5-membered heterocycle with at least one nitrogen and no more than 2 other heteroatoms."