"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    A cephalosporin is a beta-lactam antibiotic with a 6-membered dihydrothiazine ring fused to the beta-lactam ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure of cephalosporins
    # The core structure includes a beta-lactam ring fused to a 6-membered dihydrothiazine ring
    # More flexible pattern that allows for various substituents
    cephalosporin_core = Chem.MolFromSmarts("[C@H]12SCC(=C(N1C(=O)[C@@H]2*)*)*")
    
    # Check if the molecule contains the core structure
    if not mol.HasSubstructMatch(cephalosporin_core):
        return False, "No cephalosporin core structure found"

    # Check for the presence of a carboxyl group (or ester) at position 4 of the dihydrothiazine ring
    # This pattern is more flexible and allows for esterified carboxyl groups
    carboxyl_pattern = Chem.MolFromSmarts("[C@H]12SCC(=C(N1C(=O)[C@@H]2*)*)[C;$(C(=O)[OH,O])]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (or ester) at position 4 of the dihydrothiazine ring"

    return True, "Contains the core structure of a cephalosporin with a carboxyl group (or ester) at position 4"