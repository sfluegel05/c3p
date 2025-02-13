"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine has a glycerol backbone, a phosphocholine group, and two acyl chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for glycerol backbone (3 carbon chain with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Pattern for phosphocholine group: P(-[O-])(=O)OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("P([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
        
    # Pattern for ester-bonded acyl chains: C(=O)O
    # Ensure there are exactly two acyl chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} acyl chains, need exactly 2"

    return True, "Molecule is a phosphatidylcholine"