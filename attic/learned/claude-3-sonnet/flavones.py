"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: flavones
A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavone core structure broken into parts
    # Chromone core (benzopyran-4-one)
    chromone_pattern = Chem.MolFromSmarts("O=C1CC(=O)c2ccccc2O1")
    if chromone_pattern is None:
        return None, "Invalid SMARTS pattern for chromone"
    
    # Alternative more flexible chromone pattern
    alt_chromone = Chem.MolFromSmarts("O=C1C=COc2ccccc12")
    if alt_chromone is None:
        return None, "Invalid SMARTS pattern for alternative chromone"
        
    # Complete flavone core with phenyl at position 2
    flavone_core = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    if flavone_core is None:
        return None, "Invalid SMARTS pattern for flavone core"

    # Check for core structure matches
    has_core = False
    if mol.HasSubstructMatch(flavone_core):
        has_core = True
    elif mol.HasSubstructMatch(chromone_pattern) or mol.HasSubstructMatch(alt_chromone):
        # Check for attached phenyl group if we found chromone
        phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
        if phenyl_pattern is not None and mol.HasSubstructMatch(phenyl_pattern):
            has_core = True
    
    if not has_core:
        return False, "Missing flavone core structure (2-phenylchromen-4-one skeleton)"

    # Verify basic requirements
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 2:
        return False, "Must contain at least two rings"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for flavone structure"
    if o_count < 2:
        return False, f"Too few oxygens ({o_count}) for flavone structure"

    # Look for common substitutions
    substitutions = []
    
    # Common substitution patterns
    patterns = {
        "hydroxy": "cO[H]",
        "methoxy": "cOC",
        "glycoside": "[OH1,OH0][CH1]([OH1,OH0])[CH1]([OH1,OH0])[CH1,CH2]O",
        "prenyl": "CC(C)=CC"
    }
    
    for name, pattern in patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            substitutions.append(name)

    # Build response message
    msg = "Contains flavone core structure (2-phenylchromen-4-one)"
    if substitutions:
        msg += f" with {', '.join(substitutions)} substitutions"
    
    return True, msg