"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: CHEBI:6051 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes are C20 polyunsaturated fatty acids with four double bonds,
    three of which are conjugated.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons (should be around 20, allow some variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 24):
        return False, f"Carbon count {c_count} outside acceptable range for leukotriene (18-24)"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX2-]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Look for conjugated double bond systems with various possible geometries
    conjugated_patterns = [
        # Basic conjugated triene patterns with different possible geometries
        "C=CC=CC=CC",
        "C=CC=CC=CO",
        "C=CC=CC=CS",
        # Allow for different substitution patterns
        "[C,O,N,S]C=CC=CC=C[C,O,N,S]",
        # Patterns for specific leukotriene series
        "CC=CC=CC=CCO", # Common in B series
        "CC=CC=CC=CCS", # Common in C/D/E series
    ]
    
    has_conjugated = False
    for pattern in conjugated_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            has_conjugated = True
            break
    
    if not has_conjugated:
        return False, "Missing required conjugated double bond system"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_count < 4:
        return False, f"Insufficient double bonds ({double_bond_count}), need at least 4"

    # Check for characteristic features of different leukotriene series
    features = []
    
    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count > 0:
        features.append(f"{hydroxyl_count} hydroxyl groups")

    # Check for cysteinyl group (present in LTC4, LTD4, LTE4)
    cysteinyl_pattern = Chem.MolFromSmarts("SCC([NH2,NH])")
    if mol.HasSubstructMatch(cysteinyl_pattern):
        features.append("cysteinyl group")

    # Check for epoxide group (present in LTA4)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
        features.append("epoxide group")

    # Additional check for long carbon chain
    chain_pattern = Chem.MolFromSmarts("CCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing required carbon chain"

    # Verify oxygen content
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen content"

    feature_str = ", ".join(features) if features else "basic structure"
    return True, f"Leukotriene with conjugated system, {double_bond_count} double bonds, and {feature_str}"