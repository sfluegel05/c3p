"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol consists of a glycerol backbone in which one hydroxy group
    is esterified with a fatty acid, and the other two
    hydroxy groups can be either free or substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern with two hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two hydroxyl groups found"
    
    # Look for ester linkage
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Verify the presence of a single acyl group
    # There should be a carbon chain attached to the ester linkage
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)OC(C)C")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    if not fatty_acid_matches:
        return False, "Missing acyl group attached to the ester linkage"

    return True, "Contains glycerol backbone with one acyl group and two hydroxyl groups"

# Test example
# smiles = "O(C(=O)CCCCCCC/C=C\CCCCCCCC)CC(O)CO" # Example of a monoacylglycerol
# print(is_monoacylglycerol(smiles))