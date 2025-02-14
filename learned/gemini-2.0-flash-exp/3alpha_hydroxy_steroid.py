"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    Uses a SMARTS-based approach for ring detection and alpha configuration,
    targeting a specific pattern for position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for steroid core (specific tetracyclic ring system)
    # This SMARTS pattern aims to match the typical steroid core structure
    # with 3 six-membered rings and one 5-membered ring.
    steroid_core_pattern = Chem.MolFromSmarts("[C]12[C][C]([C])[C]([C])([C])3[C]([C])([C])[C]4[C]([C])([C])1[C]2[C]34")
    if steroid_core_pattern is None:
        return None, "Invalid steroid core SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_core_pattern):
          return False, "Not a steroid core structure"
      
    # 2. Check for the 3-alpha-hydroxy group using a more specific pattern, including chirality
    # This pattern assumes that the 3 carbon will have an H, a carbon, and an OH group.
    # The alpha configuration is captured by using [C@H].
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@H]1[C]([H])([O])[C]([C])([C])2[C]3[C]1([C])[C]([C])([C])([C])2[C]3")

    if alpha_hydroxy_pattern is None:
        return None, "Invalid alpha hydroxy SMARTS pattern"

    alpha_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    if len(alpha_matches) == 0:
        
        #Check for beta hydroxy group for more descriptive negative outcome
        beta_hydroxy_pattern = Chem.MolFromSmarts("[C@@H]1[C]([H])([O])[C]([C])([C])2[C]3[C]1([C])[C]([C])([C])([C])2[C]3")
        
        if beta_hydroxy_pattern is None:
           return None, "Invalid beta hydroxy SMARTS pattern"

        beta_matches = mol.GetSubstructMatches(beta_hydroxy_pattern)

        if len(beta_matches) > 0:
            return False, "Hydroxyl group at position 3 is in beta-configuration"

        return False, "No 3-alpha hydroxyl group found"

    return True, "Molecule is a 3alpha-hydroxy steroid"