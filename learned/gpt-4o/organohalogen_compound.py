"""
Classifies: CHEBI:17792 organohalogen compound
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for carbon-halogen bond
    carbon_halogen_pattern = Chem.MolFromSmarts("[CX4H0;!$(C#[!#6])]~[F,Cl,Br,I]")
    if mol.HasSubstructMatch(carbon_halogen_pattern):
        return True, "Contains carbon-halogen bond(s)"
    else:
        return False, "No carbon-halogen bond found"