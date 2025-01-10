"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Classifies a molecule as a phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated phosphate group pattern: Specifically checks for a phosphatidyl linkage
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(OC[C@H](O)COP)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphatidyl linkage found"
        
    # Updated serine pattern: accommodates stereochemistry and ester linkage
    serine_pattern = Chem.MolFromSmarts("C[C@H](N)C(=O)O[P]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine group pattern found esterified to the phosphate"

    # Glycerol backbone attached to phosphate
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No correct glycerol backbone pattern found"

    # Two ester linkages for fatty acid chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups; found {len(ester_matches)}, need at least 2"

    return True, "Molecule meets the criteria for phosphatidyl-L-serine"