"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is characterized by a glycerol backbone where one hydroxy group,
    commonly but not necessarily primary, is esterified with phosphoric acid,
    and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Phosphate group detection: P connected to three oxygens, one of which is doubly bonded
    phosphate_pat = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pat):
        return False, "Phosphate group incorrectly configured"

    # Ester group linked to glycerol with acyl chain: C(=O)OC
    ester_pat = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pat)
    
    # Look for exactly two ester linkages indicating acyl chains
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 2"

    # Glycerol backbone: checking for simple 3-carbon chain with suitable connectivity
    glycerol_pat = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pat):
        return False, "Missing glycerol backbone in expected configuration"

    return True, "Structure matches phosphatidic acid with glycerol backbone and appropriate esterifications"