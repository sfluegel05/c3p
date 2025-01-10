"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has a single fatty acid chain attached to the primary hydroxyl group of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern with two free hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two free hydroxyl groups found"
    
    # Look for a single ester linkage [C(=O)O] linked to the primary oxygen on glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Check if the ester linkage is attached to a long carbon chain (at least 3 carbon atoms)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCC([CH2])[CH2]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Ester linkage does not lead to a sufficiently long carbon chain"

    return True, "Contains glycerol backbone with a single fatty acid chain esterified at the 1-position"