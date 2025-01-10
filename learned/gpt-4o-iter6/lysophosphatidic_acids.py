"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate with specific structural features.

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

    # Identify glycerol backbone linked to phosphate: 
    # Considering patterns like [O][C@@H](CO)COP(=O)(O)O
    glycerol_phosphate_pattern = Chem.MolFromSmarts("O[C@H](CO)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol-phosphate linkage found"
        
    # Check for only one acyl chain, accounting for variations: OC(=O)C(C)*
    acyl_chain_pattern = Chem.MolFromSmarts("OC(=O)C(C)C")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl chains, need exactly 1"

    return True, "Matches all structural features of a lysophosphatidic acid"