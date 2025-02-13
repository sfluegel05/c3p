"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine has a glycerol backbone with acyl substituents at the 1- and
    2-hydroxy positions and a phosphoserine group at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern with correct stereochemistry and acyl groups
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](O[C@H]([C@@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)C)OC(=O)C)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with proper stereochemistry and acyl groups found"
    
    # Check for two ester bonds (acyl groups at 1 and 2 positions)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Incorrect number of acyl (ester) groups found. Expected at least 2, but got {len(ester_matches)}"

    # Overall match confirmation
    return True, "Matching structure for 3-sn-phosphatidyl-L-serine found"