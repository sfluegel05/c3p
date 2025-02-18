"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: CHEBI:16347 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone with two fatty acid esters and one phosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define flexible core structure pattern using SMARTS
    # Matches glycerol backbone with:
    # - Two oxygen-connected paths to ester groups (O...C=O)
    # - One phosphate group (O-P=O(O)(O))
    core_pattern = Chem.MolFromSmarts('C(O!@*C(=O))(O!@*C(=O))OP(=O)(O)O')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core glycerol-phosphate structure with two ester groups"
    
    # Verify exactly two ester groups in total (carboxylic esters)
    ester_pattern = Chem.MolFromSmarts('[OX2]C(=O)')  # Ester carbonyl oxygen
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    return True, "Contains glycerol backbone with two fatty acid esters and a phosphate group"