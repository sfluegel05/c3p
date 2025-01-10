"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the hopane skeleton using its SMILES representation
    # Hopane SMILES obtained from PubChem CID 92176
    hopane_smiles = "CC1(C)CC[C@]2(C)CC[C@H]3C[C@@]4(C)CC[C@H]5C[C@@](C)(CCCC5(C)C)CC[C@@]4(C)C3CC2C1"
    hopane_mol = Chem.MolFromSmiles(hopane_smiles)
    if hopane_mol is None:
        return False, "Failed to generate hopane skeleton structure"

    # Check if the input molecule contains the hopane skeleton
    if mol.HasSubstructMatch(hopane_mol):
        return True, "Contains hopane skeleton characteristic of hopanoids"
    else:
        return False, "Does not contain hopane skeleton"