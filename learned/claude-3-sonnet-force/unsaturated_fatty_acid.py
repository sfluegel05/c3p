"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: CHEBI:36976 unsaturated fatty acid
Any fatty acid containing at least one C=C or C#C bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for unsaturation (C=C or C#C)
    unsaturated = mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')) or mol.HasSubstructMatch(Chem.MolFromSmarts('C#C'))
    if not unsaturated:
        return False, "No unsaturation (C=C or C#C) found"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for long carbon chain (>= 8 carbons)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 8:
        return False, "Carbon chain too short (< 8 carbons)"
    
    # Check for specific patterns of unsaturation
    cis_pattern = Chem.MolFromSmarts('/C=C/')
    trans_pattern = Chem.MolFromSmarts('C=C')
    cis_matches = mol.GetSubstructMatches(cis_pattern)
    trans_matches = mol.GetSubstructMatches(trans_pattern)
    
    unsaturation_type = "unknown"
    if cis_matches:
        unsaturation_type = "cis"
    if trans_matches:
        unsaturation_type = "trans" if not cis_matches else "cis and trans"
    
    # Check for additional functional groups
    hydroxy_pattern = Chem.MolFromSmarts('O')
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    
    # Construct reason string
    reason = f"Contains {unsaturation_type} unsaturation, carboxylic acid group, and long carbon chain"
    if has_hydroxy:
        reason += ", and hydroxy group(s)"
    
    return True, reason