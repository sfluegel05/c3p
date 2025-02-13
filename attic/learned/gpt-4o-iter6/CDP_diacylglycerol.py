"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CHEBI: 
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol has a glycerol backbone with two acyl chains at positions 1 and 2, and a cytidine diphosphate group.

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

    # Look for glycerol backbone pattern (C[C@H](O)C(O))
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](O)C(O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for cytidine diphosphate group
    cdp_pattern = Chem.MolFromSmarts("[nH]1ccc(nc1=O)P(OP(=O)(O)O)(=O)OC[C@H]")
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "No cytidine diphosphate group found"

    # Look for ester linkages indicating acyl chains (-C(=O)OC-)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl chains, need at least 2"

    # Acyl chains should be long, let's check the carbon count on them
    acyl_chain_count = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and len(atom.GetBonds()) == 3 # indicative of C-C bonds in hydrocarbon chains
    )
    
    if acyl_chain_count < 20:
        return False, "Acyl chains too short"

    return True, "Contains glycerol backbone with cytidine diphosphate and two acyl chains"