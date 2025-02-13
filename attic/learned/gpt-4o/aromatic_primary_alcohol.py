"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Classifies a molecule as an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has the hydroxy group attached to a primary carbon (sp3)
    which is itself bonded directly to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Correcting the SMARTS pattern for aromatic primary alcohol
    # [CH2] is an aliphatic carbon that can be bonded to the hydroxyl group and aromatic [a].
    pattern = Chem.MolFromSmarts("[CH2][OH].[c]") 

    # Check for the presence of the aromatic primary alcohol pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains hydroxyl group attached to a methylenic carbon connected to an aromatic ring"
    
    return False, "Does not contain aromatic primary alcohol structure"