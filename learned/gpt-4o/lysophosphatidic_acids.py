"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid should contain a glycerol backbone with a phosphate group and one acyl chain linked via an ester bond.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Recognize the glycerol-phosphate backbone (allow flexibility for stereochemistry with CHO groups)
    glycerol_phosphate_fg = Chem.MolFromSmarts("O[C@H](CO)COP(=O)(O)O | O[C@@H](CO)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_fg):
        return False, "Glycerol phosphate backbone not identified"

    # Look for single ester linkage off of glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)OC[C@H](CO) | C(=O)OC[C@@H](CO)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) == 0:
        return False, "No acyl chain esterified to the glycerol backbone"
    
    if len(ester_matches) > 1:
        return False, f"Multiple ester linkages detected, count: {len(ester_matches)}; may need specific locus check"

    return True, "The molecule is a lysophosphatidic acid with an expected monoglyceride phosphate ester configuration"