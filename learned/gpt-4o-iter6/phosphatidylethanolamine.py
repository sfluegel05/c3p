"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid where a phosphatidyl group is esterified to
    the hydroxy group of ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone with correct stereochemistry (chiral center at C2)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)OCCN)OC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct stereochemistry"

    # Look for a second ester linkage on the other carbon of glycerol
    second_ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(second_ester_pattern)
    if len(ester_matches) < 2:
        return False, "Missing second ester linkage, found less than 2"

    return True, "Structure matches phosphatidylethanolamine class"

# Example usage
smiles_example = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O"
result, reason = is_phosphatidylethanolamine(smiles_example)
print(result, reason)