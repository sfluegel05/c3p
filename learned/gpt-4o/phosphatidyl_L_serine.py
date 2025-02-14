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

    # Look for chiral phosphate linkage with a serine attachment
    serine_pattern = Chem.MolFromSmarts("P(=O)(O[C@@H](N)C(=O)O)O")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No phosphoserine linkage found"
    
    # Look for two fatty acid ester chains - accounts for various carbon chain lengths and configurations
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester links, need at least 2 for fatty acid chains"
    
    # Check for glycerol backbone structure: A chain of 3 oxygen-connected carbons
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone detected"
    
    return True, "Contains glycerol backbone, two fatty acid chains, and phosphatidyl group esterified to serine"