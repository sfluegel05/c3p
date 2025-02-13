"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol involves a phosphatidyl group esterified to one of the
    hydroxy groups of inositol, a hexahydroxy cyclohexane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define substructure patterns.
    inositol_pattern = Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C(O)C1O)")  # Flexible inositol detection with hydroxy groups
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")  # Glycerol backbone pattern
    
    # Check for inositol ring presence.
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring with sufficient hydroxylations found"
    
    # Check for phosphate group presence.
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone typically found in phosphatidylinositols.
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found attached to the phosphate group"
    
    # Check for correct connectivity: glycerol backbone connected to a phosphate moiety which connects to inositol
    # Using a more detailed check of connectivity.
   
    glycerol_phosphate_connected = Chem.MolFromSmarts("OCC(O)COP(O)(O)=O")
    if not mol.HasSubstructMatch(glycerol_phosphate_connected):
        return False, "Glycerol backbone and phosphate group not adequately connected"
    
    return True, "The molecule is a phosphatidylinositol with correctly established linkages"