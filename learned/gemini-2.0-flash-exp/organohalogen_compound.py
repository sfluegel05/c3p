"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:25386 organohalogen compound
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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for a carbon-halogen bond.
    # This pattern looks for a halogen connected to a carbon atom by any bond,
    carbon_halogen_pattern = Chem.MolFromSmarts("[C]~[F,Cl,Br,I]")
    iodyl_pattern = Chem.MolFromSmarts("[C]~[I](=[O])=[O]")
    iodine_pattern = Chem.MolFromSmarts("[C]~[I]")
    
    # Check if the molecule contains the pattern
    if mol.HasSubstructMatch(carbon_halogen_pattern) or mol.HasSubstructMatch(iodyl_pattern) or mol.HasSubstructMatch(iodine_pattern):
        return True, "Contains at least one carbon-halogen bond"

    return False, "No carbon-halogen bond found"