"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:XXXXX hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid with a hopane skeleton (pentacyclic structure).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the hopane core (pentacyclic structure) without stereochemistry
    hopane_core_smarts = (
        "[C]12[C][C]3[C]([C]4[C]([C]5[C]([C](CC4)C(C)(C)CC5)CC3)CCC1C(C)(C)CC2)"
    )
    hopane_core = Chem.MolFromSmarts(hopane_core_smarts)
    if hopane_core is None:
        return None, None  # Invalid SMARTS pattern
    
    # Check if the molecule contains the hopane core
    if mol.HasSubstructMatch(hopane_core):
        return True, "Contains hopane core pentacyclic structure"
    
    # Check for hopene core with possible unsaturation
    hopene_core_smarts = (
        "[C]12[C][C]3[C]([C]4C=C[C]([C]5[C]([C](CC4)C(C)(C)CC5)CC3)CCC1C(C)(C)CC2)"
    )
    hopene_core = Chem.MolFromSmarts(hopene_core_smarts)
    if hopene_core and mol.HasSubstructMatch(hopene_core):
        return True, "Contains hopene core with unsaturation"
    
    # If no core matches, not a hopanoid
    return False, "No hopane core structure detected"