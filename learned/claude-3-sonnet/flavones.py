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
    # Note: The SMARTS pattern describes:
    # - A benzene ring fused to a pyran ring (benzopyran)
    # - A ketone at position 4
    # - A phenyl group at position 2
    flavone_core = Chem.MolFromSmarts("[#6]1([#6]2=[#6][#6]=[#6][#6]=[#6]2)=[#6][#6](=[O])[#6]2=[#6][#6]=[#6][#6]=[#6]2[O]1")
    
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Missing basic flavone core structure (2-phenylchromen-4-one)"
    
    # Check for minimum number of carbons and oxygens
    # Flavones should have at least 15 carbons and 2 oxygens in basic form
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for flavone structure"
    if o_count < 2:
        return False, f"Too few oxygens ({o_count}) for flavone structure"
    
    # Check for aromaticity
    # Both rings (A and B) should be aromatic in flavones
    ring_info = mol.GetRingInfo()
    aromatic_rings = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
            
    if aromatic_rings < 2:
        return False, "Missing required aromatic rings"
    
    # Additional check for the ketone group at position 4
    ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#6]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group at position 4"
    
    # Look for common substitution patterns (optional, but typical for flavones)
    # Hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[#6]-[OX2H1]")
    # Methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("[#6]-[OX2]-[CH3]")
    
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    
    substitution_info = []
    if has_hydroxy:
        substitution_info.append("hydroxy")
    if has_methoxy:
        substitution_info.append("methoxy")
    
    substitution_str = " and ".join(substitution_info)
    if substitution_str:
        substitution_str = f" with {substitution_str} substitutions"
    
    return True, f"Contains flavone core structure (2-phenylchromen-4-one){substitution_str}"