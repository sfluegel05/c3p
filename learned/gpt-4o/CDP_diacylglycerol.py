"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    CDP-diacylglycerol is a glycerol molecule with acyl groups (typically fatty acyl)
    at the 1- and 2-positions, and a cytidine diphosphate (CDP) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with 1,2-diacyl linkages
    diacyl_pattern = Chem.MolFromSmarts("O[C@H](COC(=O)[#6])[C@H](COC(=O)[#6])")
    if not mol.HasSubstructMatch(diacyl_pattern):
        return False, "No 1,2-diacylglycerol backbone detected"
    
    # Identify the CDP group (cytidine with diphosphate)
    cdp_pattern = Chem.MolFromSmarts("N1C=CC(=NC1=O)N")
    diphosphate_pattern = Chem.MolFromSmarts("P(=O)(O)OP(=O)(O)O")
    
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "No cytidine component"
        
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate group"
    
    # Ensure two acyl chains are present
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 2:
        return False, f"Expected 2 acyl chains but found {len(acyl_matches)}"
    
    return True, "Molecule is a CDP-diacylglycerol with the necessary components"