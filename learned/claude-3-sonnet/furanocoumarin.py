"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:47835 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic requirements
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:  # Must have at least 3 rings (benzene + pyrone + furan)
        return False, "Insufficient number of rings"

    # Count oxygens (need at least 3: one in pyrone ring, one in ether linkage, one in furan)
    o_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if o_count < 3:
        return False, "Insufficient number of oxygen atoms"

    # Basic coumarin patterns (more generalized)
    coumarin_patterns = [
        # Basic coumarin core with flexible bond types
        "[#8]-1-[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]~2-[#6](=[#8])-[#6]-1",  # Basic pattern
        "[#8]-1-[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]~2-[#6](=[#8])-[#6]=1",   # Alternative
        "[#8]-1-c~2c(C(=O)C=1)cccc2",  # Simpler pattern
        "O=C1Oc2ccccc2C1",             # Allow for different saturation
    ]
    
    has_coumarin = False
    for pat in coumarin_patterns:
        pattern = Chem.MolFromSmarts(pat)
        if pattern and mol.HasSubstructMatch(pattern):
            has_coumarin = True
            break
            
    if not has_coumarin:
        return False, "No coumarin core found"

    # Furan patterns (both saturated and unsaturated)
    furan_patterns = [
        "[#8]1[#6]~[#6]~[#6]~[#6]1",    # Basic furan/dihydrofuran
        "o1cccc1",                        # Aromatic furan
        "O1CCCC1",                        # Saturated furan
        "[#8]1-[#6]=[#6]-[#6]=[#6]1",    # Explicit double bonds
    ]
    
    has_furan = False
    for pat in furan_patterns:
        pattern = Chem.MolFromSmarts(pat)
        if pattern and mol.HasSubstructMatch(pattern):
            has_furan = True
            break
            
    if not has_furan:
        return False, "No furan ring found"

    # Check for fusion between coumarin and furan
    # These patterns look for the characteristic fusion points
    fusion_patterns = [
        # Linear fusion patterns
        "[#8]1[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~3~[#6](=[#8])~[#8]~[#6]~[#6]~3~[#6]~2",
        # Angular fusion patterns
        "[#8]1[#6]~[#6]~2~[#6]~3~[#6](=[#8])~[#8]~[#6]~[#6]~3~[#6](~[#6]~2)~[#6]~1",
        # Simplified fusion patterns
        "O1CCc2c1cc1ccc(=O)oc21",
        "O1CCc2c3c(cc12)ccc(=O)o3",
    ]
    
    has_fusion = False
    fusion_type = []
    
    for i, pat in enumerate(fusion_patterns):
        pattern = Chem.MolFromSmarts(pat)
        if pattern and mol.HasSubstructMatch(pattern):
            has_fusion = True
            if i < 2:
                fusion_type.append("linear")
            else:
                fusion_type.append("angular")
                
    if not has_fusion:
        return False, "No proper fusion between furan and coumarin rings"

    fusion_desc = f"({', '.join(set(fusion_type))} fusion)" if fusion_type else "(fusion type undetermined)"
    return True, f"Confirmed furanocoumarin structure {fusion_desc}"