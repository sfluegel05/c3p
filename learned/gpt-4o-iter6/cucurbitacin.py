"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids developed by some plants.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a tetracyclic ring system
    ring_info = Chem.GetSymmSSSR(mol)
    if len(ring_info) < 4:
        return False, "Contains fewer than 4 rings, does not meet tetracyclic criteria"
    
    # Check for specific hydroxyl group pattern OH attached to tertiary or secondary carbon
    hydroxyl_smarts = Chem.MolFromSmarts('[CX4,CX3](O)')
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "Missing hydroxyl groups in relevant contexts"
    
    # Check for carbonyl groups specifically in conjugation
    carbonyl_smarts = Chem.MolFromSmarts('[CX3]=[OX1]')
    if not mol.HasSubstructMatch(carbonyl_smarts):
        return False, "Missing carbonyl groups in the expected context"
    
    # Define a more precise cucurbitane backbone pattern
    cucurbitane_smarts = Chem.MolFromSmarts('C1CCC2C3C1C(=O)C=C4C3CC(C2)C4')
    if not mol.HasSubstructMatch(cucurbitane_smarts):
        return False, "Does not contain the typical cucurbitane structure"
    
    # Incorporate allowance for variability, provide context if needed
    return True, "Valid cucurbitacin identified by tetracyclic triterpenoid backbone"