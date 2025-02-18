"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:27102 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic benzoquinone core pattern
    # More flexible pattern that captures various substitution patterns
    benzoquinone_pattern = Chem.MolFromSmarts("[OX1]=C1C=C(C(=O)C(=C1))") 
    if benzoquinone_pattern is None:
        return False, "Invalid SMARTS pattern for benzoquinone"
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No benzoquinone core found"

    # Look for at least one methoxy group
    methoxy_pattern = Chem.MolFromSmarts("COc1c(OC|O)c(=O)c(C)c(=O)c1")
    if methoxy_pattern is None:
        return False, "Invalid SMARTS pattern for methoxy"
    if not mol.HasSubstructMatch(methoxy_pattern):
        methoxy_pattern2 = Chem.MolFromSmarts("COc1c(=O)c(C)c(OC|O)c(=O)c1")
        if methoxy_pattern2 is None or not mol.HasSubstructMatch(methoxy_pattern2):
            return False, "Missing required methoxy substitution"

    # Count oxygens (should have at least 4 for the basic structure)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Insufficient oxygen atoms for ubiquinone structure"

    # Look for characteristic side chains
    # This includes various prenyl unit patterns and longer chains
    side_chain_patterns = [
        "CC(C)=CCC",  # Basic prenyl unit
        "CC(C)CCCC",  # Reduced prenyl unit
        "CCC=C(C)C",  # Alternative prenyl arrangement
        "CCCCCCCC"    # Long alkyl chain
    ]
    
    has_side_chain = False
    for pattern in side_chain_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None and mol.HasSubstructMatch(patt):
            has_side_chain = True
            break
    
    if not has_side_chain:
        return False, "Missing characteristic side chain"

    # Additional checks for typical ubiquinone features
    
    # Check for proper ring size (should be 6-membered)
    ring_size_pattern = Chem.MolFromSmarts("[r3,r4,r5,r7,r8]")
    if ring_size_pattern is not None and mol.HasSubstructMatch(ring_size_pattern):
        return False, "Contains incorrect ring size"

    # Count carbons (should have at least 10 for minimal ubiquinone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Insufficient carbon atoms for ubiquinone structure"

    # Check for characteristic quinone oxygens pattern
    quinone_o_pattern = Chem.MolFromSmarts("O=C1[#6]~[#6]~[#6]C(=O)")
    if quinone_o_pattern is None or not mol.HasSubstructMatch(quinone_o_pattern):
        return False, "Missing proper quinone oxygen arrangement"

    return True, "Contains required ubiquinone structural features: benzoquinone core, methoxy group(s), and characteristic side chain"