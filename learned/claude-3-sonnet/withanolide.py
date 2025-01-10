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
        str: str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - withanolides typically have around 28 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (24 <= c_count <= 35):  # Wider range to accommodate derivatives
        return False, f"Carbon count {c_count} not typical for withanolide structure"

    # More flexible steroid core pattern that accounts for variations
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]4~[#6]3~[#6]2~[#6]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for characteristic δ-lactone ring in side chain
    # This pattern looks for the specific lactone arrangement found in withanolides
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1OC[C@@H]([C@H]1C)C"), # Common withanolide lactone
        Chem.MolFromSmarts("O=C1OC[CH][CH]1"), # Simpler lactone pattern
        Chem.MolFromSmarts("O=C1OC(C)=C(C)C1"), # Alternative lactone pattern
        Chem.MolFromSmarts("O=C1OC(C)(C)C=C1") # Another variation
    ]
    
    lactone_found = False
    for pattern in lactone_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            lactone_found = True
            break
    
    if not lactone_found:
        return False, "No characteristic withanolide lactone ring found"

    # Check for oxygen count (including lactone oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if not (2 <= o_count <= 15):  # Wider range to accommodate derivatives and glycosides
        return False, f"Oxygen count {o_count} not typical for withanolides"

    # Check for reasonable molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 900):  # Wider range to accommodate derivatives
        return False, f"Molecular weight {mol_wt} outside typical range for withanolides"

    # Check for ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:
        return False, f"Ring count {ring_count} too low for withanolide structure"

    # Additional structural features common in withanolides
    # Look for characteristic oxygen substitutions
    hydroxy_pattern = Chem.MolFromSmarts("[#6]~[OX2H1]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Missing typical oxygen substitutions"

    # If all checks pass, it's likely a withanolide
    return True, "Matches withanolide structure with steroid core, lactone ring, and characteristic substitution pattern"