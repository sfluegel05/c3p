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
    
    # Get basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for glycosylated/esterified carotenoids
    is_glycosylated = mol.HasSubstructMatch(Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1"))
    is_esterified = mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)OC"))
    
    # Adjust carbon count limits based on modifications
    max_carbons = 50
    if is_glycosylated:
        max_carbons = 80  # Allow for glycosylated forms
    if is_esterified:
        max_carbons = 60  # Allow for esterified forms
        
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) for a carotenoid"
    if c_count > max_carbons:
        return False, f"Too many carbons ({c_count}) for a carotenoid"
    
    # Check for carotenoid backbone patterns
    backbone_patterns = [
        "C(/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C))",  # All-trans polyene
        "C(/C=C/C=C(/C)C=C/C=C(/C))",                # Shorter polyene
        "CC(C)=CCC/C(C)=C/CC/C(C)=C/CC"              # Carotenoid precursor
    ]
    
    has_backbone = False
    for pattern in backbone_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_backbone = True
            break
            
    if not has_backbone:
        return False, "No characteristic carotenoid backbone found"
    
    # Count conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    num_conjugated = len(mol.GetSubstructMatches(conjugated_pattern))
    
    # Look for common end groups with more specific patterns
    end_groups = [
        ("beta", "C1=C(C)C(C)(C)CCC1"),           # Beta-type end
        ("epsilon", "C1C(C)=C(C)CC1(C)C"),        # Epsilon-type end
        ("gamma", "C1=C(C)CCCC1(C)C"),            # Gamma-type end
        ("keto", "C(=O)C=CC=C"),                  # Keto end group
        ("acyclic", "C(C)(C)C=CC=C"),             # Acyclic end
        ("cyclopentyl", "C1CCCC1(C)C"),           # Cyclopentyl end
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
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):
        modifications.append("hydroxylated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)C")):
        modifications.append("keto")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        modifications.append("epoxidated")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O")):
        modifications.append("carboxylated")
    
    # Count methyl branches (characteristic of carotenoids)
    methyl_pattern = Chem.MolFromSmarts("C-C(C)-C")
    num_methyls = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Final classification logic
    is_likely_carotenoid = False
    reason = ""
    
    if found_ends and num_conjugated >= 4:
        is_likely_carotenoid = True
        mod_str = f" ({', '.join(modifications)})" if modifications else ""
        end_str = ", ".join(found_ends)
        reason = f"Carotenoid with {end_str} end group(s){mod_str}"
    elif num_conjugated >= 7 and num_methyls >= 4:
        is_likely_carotenoid = True
        reason = "Carotenoid based on conjugation pattern and methyl branching"
    elif c_count >= 30 and num_conjugated >= 6 and has_backbone:
        is_likely_carotenoid = True
        reason = "Carotenoid based on backbone structure"
    else:
        reason = "Does not match carotenoid structural requirements"
        
    return is_likely_carotenoid, reason