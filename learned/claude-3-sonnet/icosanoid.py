"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:24913 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count oxygens (should have multiple from oxidation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Too few oxygen atoms for an icosanoid"

    # Look for key structural features
    features_found = []
    
    # Carbon chain patterns characteristic of icosanoids
    chain_patterns = [
        'C~C~C~C~C~C~C~C~C~C~C~C', # Basic long chain
        'C~C=C~C=C~C~C', # Conjugated double bonds
        'C1CCCC1C~C~C~C', # Prostaglandin-like
    ]
    
    for pattern in chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            features_found.append("characteristic carbon chain")
            break

    if not features_found:
        return False, "No characteristic icosanoid carbon chain found"

    # Oxidation and functional group patterns
    patterns = {
        'C=C': "double bonds",
        '[OH]': "hydroxyl groups",
        'C1CCCC1': "cyclopentane ring",
        'C1OC1': "epoxide group",
        'C(=O)C': "ketone group",
        'C(=O)[OH]': "carboxylic acid",
        'C(=O)O[CH2,CH3]': "ester group",
        '[OO]': "peroxide group",
        'C(=O)N': "amide group",
    }
    
    for smart, feature in patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smart)):
            features_found.append(feature)

    # Need multiple oxidation features to be an icosanoid
    oxidation_features = ["hydroxyl groups", "epoxide group", "ketone group", 
                         "carboxylic acid", "peroxide group"]
    oxidation_count = sum(1 for f in oxidation_features if f in features_found)
    
    if oxidation_count < 1:
        return False, "Lacks typical icosanoid oxidation features"

    # Look for specific icosanoid substructures
    icosanoid_patterns = [
        'C(=O)CCC=CCC1C(O)CC(O)C1', # Prostaglandin core
        'C=CC=CC=CC(O)', # Characteristic hydroxylated polyene
        'CC=CC=CC=CC(OO)', # Hydroperoxy pattern
        'C1OC1CC=CC=C', # Epoxy-fatty acid pattern
    ]
    
    for pattern in icosanoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            features_found.append("characteristic icosanoid substructure")
            break

    # If we have enough characteristic features, classify as icosanoid
    if len(set(features_found)) >= 3:  # Need at least 3 different features
        features_str = ", ".join(set(features_found))
        return True, f"Contains characteristic icosanoid features: {features_str}"
    
    return False, "Insufficient characteristic icosanoid features"