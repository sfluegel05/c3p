"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    CDP-diacylglycerol is defined as a glycerol molecule with acyl groups (typically fatty acyl)
    at the 1- and 2-positions, and a cytidine diphosphate (CDP) group at the 3-position.

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
    
    # Identify the 1,2-diacylglycerol backbone
    diacylglycerol_pattern = Chem.MolFromSmarts("O[C@H](COP(=O)(O)O[P](=O)(O)O)COC(=O)")
    if not mol.HasSubstructMatch(diacylglycerol_pattern):
        return False, "No 1,2-diacylglycerol backbone with CDP group at 3-position"
    
    # Identify the cytidine component
    cytidine_pattern = Chem.MolFromSmarts("N1C=CC(=NC1=O)N")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine component present in structure"
    
    # Check for two acyl chains (esters typically imply this)
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_chain_matches) < 2:
        return False, f"Less than two acyl (ester) chains found, got {len(acyl_chain_matches)}"

    return True, "Contains characteristic CDP-diacylglycerol backbone with cytidine and acyl groups"