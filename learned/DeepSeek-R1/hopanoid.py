"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:XXXXX hopanoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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
    
    # Define SMARTS pattern for the hopane core (pentacyclic structure)
    hopane_core = Chem.MolFromSmarts(
        "[C@]12[C@@H]3[C@@H]([C@@H]4[C@]([C@@H]([C@]5([C@@H](CC4)C(C)(C)CC5)C)CC3)C)CC[C@H]1C(C)(C)CC2"
    )
    if hopane_core is None:
        return None, None  # Invalid SMARTS pattern
    
    # Check if the molecule contains the hopane core
    if mol.HasSubstructMatch(hopane_core):
        return True, "Contains hopane core pentacyclic structure"
    
    # Additional check for variants with unsaturation in the core
    hopene_core = Chem.MolFromSmarts(
        "[C@]12[C@@H]3[C@@H]([C@@H]4C=C[C@]([C@@H]([C@]5([C@@H](CC4)C(C)(C)CC5)C)CC3)C)CC[C@H]1C(C)(C)CC2"
    )
    if hopene_core and mol.HasSubstructMatch(hopene_core):
        return True, "Contains hopene core with unsaturation"
    
    # If no core matches, not a hopanoid
    return False, "No hopane core structure detected"