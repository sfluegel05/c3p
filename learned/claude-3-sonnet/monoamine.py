"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine has one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a1aaaaa1")  # 6-membered aromatic ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"

    # Look for primary, secondary, or tertiary amines
    # [#7X3] matches any trivalent nitrogen (typical amine)
    # The !$(N=[*]) part excludes imines and similar
    amine_pattern = Chem.MolFromSmarts("[#7X3;!$(N=[*])]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    
    if not amine_matches:
        return False, "No amine group found"
    
    # Common monoamine patterns with 2-carbon linker between amine and aromatic ring
    # Pattern 1: -CH2-CH2-N (ethylamine)
    # Pattern 2: -CH(OH)-CH2-N (β-hydroxylethylamine, like in adrenaline)
    # Pattern 3: -CH2-CH(NH2)- (like in some amino acids)
    monoamine_patterns = [
        Chem.MolFromSmarts("a-[CH2][CH2][NX3]"), # simple ethylamine linker
        Chem.MolFromSmarts("a-[CH2][CH]([OH])[NX3]"), # beta-hydroxyl pattern
        Chem.MolFromSmarts("a-[CH]([OH])[CH2][NX3]"), # alternate hydroxyl pattern
        Chem.MolFromSmarts("a-[CH2][CH]([NH2])-"), # amino acid like pattern
    ]
    
    found_valid_pattern = False
    for pattern in monoamine_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_valid_pattern = True
            break
            
    if not found_valid_pattern:
        return False, "No valid two-carbon chain connecting amine to aromatic ring"
    
    # Additional check for common monoamine features
    # Look for hydroxyl groups on aromatic ring (common in natural monoamines)
    aromatic_oh_pattern = Chem.MolFromSmarts("aO[H]")
    has_aromatic_oh = mol.HasSubstructMatch(aromatic_oh_pattern)
    
    # Count number of amine groups to ensure we don't have too many
    amine_count = len(mol.GetSubstructMatches(amine_pattern))
    if amine_count > 3:  # allowing up to 3 for some derivatives
        return False, "Too many amine groups"
    
    # Success message with additional detail
    detail = "monoamine with"
    if has_aromatic_oh:
        detail += " hydroxylated aromatic ring"
    else:
        detail += " aromatic ring"
    
    return True, f"Confirmed {detail} and two-carbon amine linker"