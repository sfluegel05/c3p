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

    # Basic flavone core structure (2-phenylchromen-4-one)
    # More flexible pattern that accounts for aromatic bonds and substitutions
    flavone_core = Chem.MolFromSmarts("[#6]1=[#6]-c2c([#6](=[O])-[#6]1)c([*,H])c([*,H])c([*,H])c2[*,H]")
    
    # Alternative pattern for the core
    alt_core = Chem.MolFromSmarts("O=C1C=C(Oc2ccccc12)c1ccc([*,H])cc1")
    
    # Pattern for checking the essential oxygen arrangement
    oxygen_pattern = Chem.MolFromSmarts("O=C1[#6]~[#6]~[#6]Oc2ccccc12")

    if not (flavone_core and alt_core and oxygen_pattern):
        return None, "Invalid SMARTS patterns"

    # Check for core structure matches
    has_core = False
    if mol.HasSubstructMatch(flavone_core) or mol.HasSubstructMatch(alt_core):
        if mol.HasSubstructMatch(oxygen_pattern):
            has_core = True

    if not has_core:
        return False, "Missing flavone core structure (2-phenylchromen-4-one skeleton)"

    # Verify ring structure
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 2:
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
    
    patterns = {
        "hydroxy": "[OH1]",
        "methoxy": "CO[c]",
        "glycoside": "[OH1,OH0][CH1]([OH1,OH0])[CH1]([OH1,OH0])[CH1,CH2]O",
        "prenyl": "CC(C)=CC[c]",
        "sulfate": "OS(=O)(=O)[OH]"
    }
    
    for name, pattern in patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            substitutions.append(name)

    # Additional check for aromaticity
    aromatic_rings = sum(1 for ring in rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if aromatic_rings < 2:
        return False, "Must contain at least two aromatic rings"

    # Build response message
    msg = "Contains flavone core structure (2-phenylchromen-4-one)"
    if substitutions:
        msg += f" with {', '.join(substitutions)} substitutions"
    
    return True, msg