"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: CHEBI:36584 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is a carnitine molecule with an acyl group attached to the
    oxygen atom, and the carnitine component must have the L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carnitine backbone pattern ([N+]([C])([C])[C@H](CC([O-])=O))
    carnitine_pattern = Chem.MolFromSmarts("[N+]([C])([C])[C@H](CC([O-])=O)")
    carnitine_match = mol.GetSubstructMatch(carnitine_pattern)
    if not carnitine_match:
        return False, "No L-carnitine backbone found"

    # Get the oxygen atom index of the carnitine backbone
    carnitine_oxygen_idx = carnitine_match[3]

    # Look for acyl group attached to the carnitine oxygen (-O-C(=O)-)
    acyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, expected exactly 1"

    acyl_oxygen_idx = acyl_matches[0][0]
    if acyl_oxygen_idx != carnitine_oxygen_idx:
        return False, "Acyl group not attached to the oxygen of the carnitine backbone"

    return True, "Contains L-carnitine backbone with an acyl group attached to the oxygen"