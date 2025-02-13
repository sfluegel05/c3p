"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the structure, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycerol backbone with stereochemistry [C@H](COC)(O)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COP)(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No correct glycerol backbone found"
    
    # Look for ether linkage to identify alkyl chain at position 1
    ether_linkage = Chem.MolFromSmarts("OCC")
    if not mol.HasSubstructMatch(ether_linkage):
        return False, "No ether linkage for alkyl chain found"

    # Look for ester linkage to identify acyl chain at position 2
    ester_linkage = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_linkage):
        return False, "No ester linkage for acyl chain found"
    
    # Phosphocholine group verification
    phosphocholine_group = Chem.MolFromSmarts("N(C)(C)CCOP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphocholine_group):
        return False, "No phosphocholine group found"

    # If all checks pass
    return True, "Structure matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"