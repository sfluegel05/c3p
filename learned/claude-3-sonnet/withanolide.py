"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a C28 steroid lactone with a modified side chain forming 
    a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - withanolides have C28 backbone
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (26 <= c_count <= 32):  # Allow some variation for substitutions
        return False, f"Carbon count {c_count} not typical for withanolide structure (expect ~28)"

    # More flexible steroid core patterns
    steroid_patterns = [
        # Basic 4-ring steroid system with flexible connections
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        # Alternative pattern with more flexibility
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~[#6]~1",
        # Pattern focusing on ring connectivity
        "[#6]~1~2~[#6]~[#6]~[#6]~[#6]~1~[#6]~1~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~3~[#6]~1~[#6]~2"
    ]
    
    core_found = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            core_found = True
            break
            
    if not core_found:
        return False, "No characteristic steroid core structure found"

    # Specific withanolide lactone patterns
    lactone_patterns = [
        # α,β-unsaturated δ-lactone
        "O=C1OC[C,C=][C,C=]1",
        # More specific withanolide lactone patterns
        "O=C1OC(C)=C(C)C1",
        "O=C1OCC(C)=C1C",
        # General lactone patterns with substitutions
        "O=C1OC([#6])[#6]C1([#6])[#6]",
        "O=C1OC([#6])=C([#6])C1",
        # Side chain lactone connection
        "[#6]~2~[#6]~[#6]~1~O~C(=O)~[#6]~[#6]~1~[#6]~2"
    ]
    
    lactone_found = False
    for pattern in lactone_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            lactone_found = True
            break
    
    if not lactone_found:
        return False, "No characteristic lactone ring found"

    # Check for characteristic oxygen substitutions
    oxygen_patterns = [
        "[OX2H1]-[#6]~[#6]~[#6]", # hydroxyl groups
        "[#6]-[CX3](=O)-[#6]", # ketone groups
        "[#6]~[#6](=O)~O~[#6]", # ester groups
        "[OX2H1]-[#6]-[#6]-[OX2H1]", # 1,2-diol pattern
        "O=C-C=C", # α,β-unsaturated carbonyl
    ]
    
    oxygen_features = 0
    for pattern in oxygen_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            oxygen_features += 1
    
    if oxygen_features < 2:
        return False, "Missing characteristic oxygen substitution pattern"

    # Check ring count
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if not (5 <= ring_count <= 8):
        return False, f"Ring count {ring_count} not typical for withanolide structure"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (400 <= mol_wt <= 800):
        return False, f"Molecular weight {mol_wt} outside typical range for withanolides"

    # If all checks pass
    return True, "Matches withanolide structure with characteristic steroid core, lactone ring, and substitution pattern"