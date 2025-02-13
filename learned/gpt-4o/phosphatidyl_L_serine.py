"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine consists of a glycerol backbone with two fatty acid chains,
    and a phosphatidyl group esterified to the hydroxy group of serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with phosphorus and serine esterification allowing stereochemistry variation
    structural_pattern = Chem.MolFromSmarts("OCC(O[P](=O)(O)OCC(N)C(=O)O)COP(O)(=O)O")
    
    if not mol.HasSubstructMatch(structural_pattern):
        return False, "No phosphatidyl-L-serine backbone structure detected"
    
    # Find ester groups connected with carbon chains, accounting for potential variations
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Ensure there are at least two ester linkages denoting two fatty acid chains
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} fatty acid ester links, need at least 2"
    
    return True, "Contains glycerol backbone with two fatty acid chains and phosphatidyl group esterified to serine"