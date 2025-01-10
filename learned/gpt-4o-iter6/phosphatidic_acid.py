"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid features: 
    - A glycerol backbone 
    - Two esterified fatty acids 
    - A phosphate group esterified to glycerol
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphate group linked to glycerol
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)OCCO")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group incorrectly configured"

    # Identifying two ester linkages connected to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)OC(CO)CO")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester linkages, found {len(ester_matches)}"

    # Extra validation for glycerol backbone and ester configurations
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone missing or improperly configured"

    return True, "Molecule is classified as phosphatidic acid with correct glycerol backbone and esterifications"