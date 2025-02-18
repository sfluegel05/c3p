"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine (CHEBI:74544)
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has a carnitine backbone with L-configuration and an O-acyl group attached via ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the core structure:
    # - Trimethylammonium group ([N+](C)(C)C)
    # - Chiral center (R configuration) connected to:
    #   - Ammonium group
    #   - Ester oxygen (OC(=O))
    #   - Carboxylate group (CC(=O)[O-])
    core_smarts = Chem.MolFromSmarts('[N+](C)(C)C[C@H](CC(=O)[O-])OC(=O)')
    if not mol.HasSubstructMatch(core_smarts):
        return False, "Does not match O-acyl-L-carnitine core structure"
    
    return True, "Matches O-acyl-L-carnitine core structure"