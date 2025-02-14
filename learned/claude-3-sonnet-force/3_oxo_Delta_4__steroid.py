"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:17400 3-oxo-Delta(4) steroid
A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for 3-oxo pattern
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)")
    oxo_match = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_match:
        return False, "No 3-oxo group found"
    
    # Look for Delta(4) pattern (C=C at positions 4 and 5 in a ring)
    delta4_pattern = Chem.MolFromSmarts("[CR2]=C[CR2]")
    delta4_match = mol.GetSubstructMatches(delta4_pattern)
    if not delta4_match:
        return False, "No Delta(4) alkene found"
    
    # Check if the Delta(4) alkene and 3-oxo group are part of a fused ring system
    ring_info = mol.GetRingInfo()
    is_fused = False
    for ring in ring_info.AtomRings():
        if len(ring) >= 6:
            ring_atoms = set(ring)
            if oxo_match[0][0] in ring_atoms and any(delta4_match[0][i] in ring_atoms for i in range(2)):
                is_fused = True
                break
    
    if not is_fused:
        return False, "Delta(4) alkene and 3-oxo group not part of a fused ring system"
    
    # Check for steroid skeleton
    steroid_pattern = Chem.MolFromSmarts("[CR4]1[CR3]@[CR3]@[CR3]@[CR3]@[CR3]@[CR3]@1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid skeleton found"
    
    return True, "Contains a 3-oxo group and a Delta(4) alkene in a fused ring system, steroid skeleton present"