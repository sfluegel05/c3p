"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: germacranolide
A sesquiterpene lactone based on germacrane skeleton
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for characteristic α-methylene-γ-lactone patterns
    lactone_patterns = [
        # α-methylene-γ-lactone patterns with various connections to 10-membered ring
        Chem.MolFromSmarts("O=C1OC2CC[C@@H]1*2"),  # Basic fused pattern
        Chem.MolFromSmarts("O=C1OC2CCC1C2"),       # Alternative connection
        Chem.MolFromSmarts("C=C1C(=O)OC2CC1*2"),   # With exocyclic methylene
        Chem.MolFromSmarts("O=C1OC([CH2])C1"),     # Simple α-methylene-γ-lactone
        # Additional patterns to catch variations
        Chem.MolFromSmarts("O=C1OC(=C)C1"),
        Chem.MolFromSmarts("O=C1OC(C=C)C1")
    ]
    
    has_lactone = False
    for pat in lactone_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_lactone = True
            break
            
    if not has_lactone:
        return False, "No characteristic α-methylene-γ-lactone found"

    # Check for 10-membered ring (germacrane skeleton)
    # More specific patterns that account for common double bonds and substitutions
    germacrane_patterns = [
        Chem.MolFromSmarts("C1CC=C(C)CCC(C)=CCC1"),  # Common germacrane pattern
        Chem.MolFromSmarts("C1CCCCCCCCC1"),          # Basic 10-membered ring
        Chem.MolFromSmarts("C1CC=CCCC=CCCC1"),       # With typical double bonds
        Chem.MolFromSmarts("C1CC(C)=CCC(C)=CCCC1")   # With methyl substitutions
    ]
    
    has_germacrane = False
    for pat in germacrane_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_germacrane = True
            break

    if not has_germacrane:
        return False, "No germacrane (10-membered ring) skeleton found"

    # Count carbons (should be ~15 for sesquiterpene core)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13:
        return False, f"Carbon count {c_count} too low for sesquiterpene lactone"
    
    # Count oxygens (should have at least 2 for lactone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for lactone structure"

    # Check for characteristic features
    features = []
    
    # Check for additional functional groups often present
    if mol.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)O")):
        features.append("acetate group")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
        features.append("hydroxyl group")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        features.append("double bonds")
    
    # Count rings (should have at least 2: 10-membered + lactone)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring count for germacranolide structure"

    # Build reason string
    reason = "Contains germacrane skeleton with fused α-methylene-γ-lactone"
    if features:
        reason += f" and {', '.join(features)}"

    return True, reason