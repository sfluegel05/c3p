"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: Flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone in which the hydrogen at position 3 of the heterocyclic ring 
    is replaced by a hydroxy group (i.e., it's a 3-hydroxyflavone).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 3-hydroxyflavone core SMARTS pattern with aromatic bonds
    flavonol_smarts = 'c1cc2c(c(c1)O)oc(=O)cc2'  # Aromatic rings with OH at position 3 and ketone at position 4

    flavonol_pattern = Chem.MolFromSmarts(flavonol_smarts)
    if flavonol_pattern is None:
        return False, "Invalid SMARTS pattern for flavonol core"

    # Check for substructure match
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "Does not contain the 3-hydroxyflavone core structure"

    return True, "Contains the 3-hydroxyflavone core structure characteristic of flavonols"