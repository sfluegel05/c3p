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
    Excludes retinoids and requires specific structural features.
    
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
    
    # Get basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for modifications that affect carbon count limits
    is_glycosylated = mol.HasSubstructMatch(Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1"))
    is_esterified = mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)OC"))
    has_phosphate = mol.HasSubstructMatch(Chem.MolFromSmarts("P(=O)(O)(O)O"))
    
    # Adjust carbon count limits based on modifications
    base_min_carbons = 30  # C30 minimum for basic carotenoids
    max_carbons = 120 if (is_glycosylated or is_esterified or has_phosphate) else 50
    
    if c_count < base_min_carbons:
        return False, f"Too few carbons ({c_count}) for a carotenoid"
    if c_count > max_carbons:
        return False, f"Too many carbons ({c_count}) for a carotenoid"
    
    # Check for retinoid patterns (to exclude)
    retinoid_patterns = [
        "CC(C)=CCC\C(C)=C\C=C\C(C)=C\C(=O)O",  # Retinoic acid
        "CC(C)=CCC\C(C)=C\C=C\C(C)=C\C=O",     # Retinal
        "CC(C)=CCC\C(C)=C\C=C\C(C)=C\CO"       # Retinol
    ]
    
    for pattern in retinoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Molecule appears to be a retinoid"
    
    # Check for carotenoid backbone patterns (more specific)
    backbone_patterns = [
        # Long conjugated polyene chain with methyl branches
        "CC(/C=C/C=C(C)/C=C/C=C(C)/C=C/C)",
        # Characteristic end-to-end pattern
        "C(/C=C/C=C(C)/C=C/C=C(C)/C=C/C=C(C)/C=C/C)",
        # Pattern with cyclic end groups
        "C1C(C)(C)CCC=C1/C=C/C=C(C)/C=C/C"
    ]
    
    has_backbone = False
    for pattern in backbone_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_backbone = True
            break
            
    if not has_backbone:
        return False, "No characteristic carotenoid backbone found"
    
    # Count conjugated double bonds (require longer conjugation)
    conjugated_pattern = Chem.MolFromSmarts("C=CC=CC=CC=C")
    num_conjugated = len(mol.GetSubstructMatches(conjugated_pattern))
    
    if num_conjugated < 2:  # Need at least two sets of 4 conjugated bonds
        return False, "Insufficient conjugation for a carotenoid"
    
    # Look for characteristic end groups
    end_groups = [
        ("beta", "C1=C(C)C(C)(C)CCC1"),
        ("epsilon", "C1C(C)=CCCC1(C)C"),
        ("gamma", "C1=C(C)CCCC1(C)C"),
        ("keto", "C(=O)C=CC=C"),
        ("cyclopentyl", "C1CCCC1(C)C")
    ]
    
    found_ends = []
    for name, smarts in end_groups:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_ends.append(name)
    
    # Check for modifications
    modifications = []
    if is_glycosylated:
        modifications.append("glycosylated")
    if is_esterified:
        modifications.append("esterified")
    if has_phosphate:
        modifications.append("phosphorylated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):
        modifications.append("hydroxylated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)C")):
        modifications.append("keto")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        modifications.append("epoxidated")
    
    # Count methyl branches (characteristic of carotenoids)
    methyl_pattern = Chem.MolFromSmarts("CC(C)C")
    num_methyls = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Final classification logic
    if found_ends and num_conjugated >= 2 and num_methyls >= 4:
        mod_str = f" ({', '.join(modifications)})" if modifications else ""
        end_str = ", ".join(found_ends)
        return True, f"Carotenoid with {end_str} end group(s){mod_str}"
    elif num_conjugated >= 3 and num_methyls >= 6 and has_backbone:
        return True, "Carotenoid based on conjugation pattern and methyl branching"
    elif c_count >= 40 and has_backbone and (is_glycosylated or has_phosphate):
        return True, f"Modified carotenoid ({', '.join(modifications)})"
    
    return False, "Does not match carotenoid structural requirements"