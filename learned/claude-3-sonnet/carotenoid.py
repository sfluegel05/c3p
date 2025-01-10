"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23042 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons - carotenoids typically have ~40 carbons but can be degraded
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Allow for apocarotenoids
        return False, f"Too few carbons ({c_count}) for a carotenoid"
    if c_count > 50:  # Allow for some modifications but not too many
        return False, f"Too many carbons ({c_count}) for a carotenoid"
    
    # Look for extended polyene chain patterns
    polyene_pattern1 = Chem.MolFromSmarts("C=CC=CC=CC=CC=C")  # Longer pattern
    polyene_pattern2 = Chem.MolFromSmarts("C=CC=CC=CC=C")     # Shorter pattern
    
    if not (mol.HasSubstructMatch(polyene_pattern1) or 
            len(mol.GetSubstructMatches(polyene_pattern2)) >= 2):
        return False, "No characteristic polyene chain found"
    
    # Count conjugated double bonds using SMARTS
    db_pattern = Chem.MolFromSmarts("C=C")
    num_double_bonds = len(mol.GetSubstructMatches(db_pattern))
    if num_double_bonds < 7:  # Minimum for most carotenoids
        return False, f"Insufficient conjugation ({num_double_bonds} double bonds)"
    
    # Look for common end groups
    end_groups = [
        ("beta", "C1C=C(C)CCC1(C)C"),          # Beta-type end
        ("epsilon", "C1C=C(C)CC1(C)C"),         # Epsilon-type end
        ("gamma", "C1C(C)=CCCC1(C)C"),         # Gamma-type end
        ("acyclic", "CC(C)=CCCC(C)(C)O"),      # Acyclic oxygenated end
        ("keto", "CC(=O)C=C"),                 # Keto end group
    ]
    
    found_ends = []
    for name, smarts in end_groups:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_ends.append(name)
    
    # Check for common modifications
    modifications = []
    
    # Hydroxyl groups (xanthophylls)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):
        modifications.append("hydroxylated")
    
    # Keto groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)C")):
        modifications.append("keto")
    
    # Epoxide groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        modifications.append("epoxidated")
    
    # Carboxylic acid groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O")):
        modifications.append("carboxylated")
    
    # Final classification
    if not found_ends:
        if c_count >= 30 and num_double_bonds >= 9:
            # Might be a carotenoid with unusual end groups
            return True, "Likely carotenoid based on carbon count and conjugation"
        return False, "No characteristic end groups found"
    
    mod_str = " (" + ", ".join(modifications) + ")" if modifications else ""
    end_str = ", ".join(found_ends)
    return True, f"Carotenoid with {end_str} end group(s){mod_str}"