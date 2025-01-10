"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:78512 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are biologically active signalling molecules made by 
    oxygenation of C20 fatty acids, excluding classic icosanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatches(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Count carbons in main chain
    carbon_chain = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No C20 carbon chain found"

    # Check for classic icosanoid patterns
    classic_patterns = [
        (Chem.MolFromSmarts("C1CC(O)C(=O)C1"), "prostaglandin-like cyclopentane ring"),
        (Chem.MolFromSmarts("C1OC1CCC1OC1"), "thromboxane-like bicycle"),
        (Chem.MolFromSmarts("[CH2][CH2][CH]1[CH2][CH2][CH]([CH2][CH2][CH2]C(=O)[OH])[CH]1"), "prostaglandin core"),
    ]
    
    for pattern, reason in classic_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, f"Classic icosanoid detected: {reason}"

    # Look for characteristic nonclassic icosanoid features
    features = []
    
    # Count oxygenated functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OH1]")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))
    epoxy_count = len(mol.GetSubstructMatches(epoxy_pattern))
    
    if hydroxy_count == 0 and epoxy_count == 0:
        return False, "No hydroxyl or epoxy groups found"
    
    # Check for conjugated double bond systems
    conjugated_patterns = [
        Chem.MolFromSmarts("C=CC=CC=C"),  # Three conjugated
        Chem.MolFromSmarts("C=CC=C"),      # Two conjugated
    ]
    
    conjugated_count = 0
    for pattern in conjugated_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            conjugated_count += len(mol.GetSubstructMatches(pattern))
    
    # Collect structural features
    if conjugated_count > 0:
        features.append(f"conjugated double bond system")
    if hydroxy_count > 0:
        features.append(f"{hydroxy_count} hydroxyl groups")
    if epoxy_count > 0:
        features.append(f"{epoxy_count} epoxy groups")

    # Additional checks for specific nonclassic icosanoid characteristics
    if hydroxy_count >= 2 and conjugated_count > 0:
        features.append("multiple hydroxylation with conjugated system")
    elif epoxy_count >= 1 and hydroxy_count >= 1:
        features.append("epoxy-hydroxy combination")
    else:
        return False, "Insufficient oxygenation pattern for nonclassic icosanoid"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low"

    # Count total oxygens (should have multiple for oxygenation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for oxygenated fatty acid"

    return True, "Nonclassic icosanoid: " + ", ".join(features)