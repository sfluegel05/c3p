"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI:17855 N-acylglycine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is a glycine molecule where the amino group is acylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycine backbone pattern (N-CH2-COOH), where N may be acylated
    glycine_pattern = Chem.MolFromSmarts("[NX3][CH2]C(=O)[OH]")
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine backbone found"

    # Look for acyl group (R-CO-) attached to the nitrogen
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) == 0:
        return False, "No acyl group found attached to the nitrogen"

    # Check if the acyl group is attached to the nitrogen of the glycine backbone
    glycine_match = mol.GetSubstructMatch(glycine_pattern)
    glycine_nitrogen_idx = glycine_match[0]  # Index of the nitrogen in the glycine backbone
    for match in acyl_matches:
        if match[1] == glycine_nitrogen_idx:  # Check if the nitrogen in the acyl group matches the glycine nitrogen
            return True, "Contains glycine backbone with acyl group attached to nitrogen"

    return False, "Acyl group not attached to the nitrogen of the glycine backbone"