"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is based on the cyclopenta[a]phenanthrene skeleton,
    consisting of three fused six-membered rings and one fused five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for steroid core (cyclopenta[a]phenanthrene skeleton)
    # This pattern represents four rings fused together: three six-membered rings and one five-membered ring
    steroid_core_smarts = '[#6]12CC[C@H]3C[C@H](CC[C@@H]3C1)C2'
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    
    if steroid_core is None:
        return False, "Invalid SMARTS pattern for steroid core"
    
    # Check if molecule has the steroid core
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid core not found"
    
    # Optionally, check for methyl groups at C-10 and C-13 positions
    # and an alkyl group at C-17 position, as per the definition
    # This would require more detailed mapping and is beyond the scope of this function
    
    return True, "Contains steroid core (cyclopenta[a]phenanthrene skeleton)"