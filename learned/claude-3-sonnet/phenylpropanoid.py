"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:35697 phenylpropanoid

A phenylpropanoid is any organic aromatic compound with a structure based on a phenylpropane skeleton.
The class includes naturally occurring phenylpropanoid esters, flavonoids, anthocyanins, coumarins and 
many small phenolic molecules as well as their semi-synthetic and synthetic analogues. Phenylpropanoids 
are also precursors of lignin.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenylpropane skeleton (C6-C3 pattern)
    phenylpropane_patterns = [
        Chem.MolFromSmarts("[c:1]1[c:2][c:3][c:4][c:5][c:6][c:7]1[C:8][C:9][C:10]"),
        Chem.MolFromSmarts("[c:1]1[c:2][c:3][c:4][c:5][c:6][c:7]1[C:8]=[C:9][C:10]")
    ]
    
    phenylpropane_match = False
    for pattern in phenylpropane_patterns:
        if mol.HasSubstructMatch(pattern):
            phenylpropane_match = True
            break
    
    if not phenylpropane_match:
        return False, "No phenylpropane skeleton found"
    
    # Check for aromatic ring
    aromatic_rings = mol.GetAromaticRings()
    if not aromatic_rings:
        return False, "No aromatic ring found"
    
    # Check for oxy substitutions (OH, OR, O-glycosides)
    oxy_pattern = Chem.MolFromSmarts("[OX2H,OX2R,OX2c]")
    oxy_matches = mol.GetSubstructMatches(oxy_pattern)
    if not oxy_matches:
        return False, "No oxy substitutions found"
    
    # Check for other common phenylpropanoid features (esters, ethers, lactones)
    feature_patterns = [
        Chem.MolFromSmarts("[OX2][CX3](=[OX1])"), # Esters
        Chem.MolFromSmarts("[OX2R0]"), # Ethers
        Chem.MolFromSmarts("[OX2r5]")  # Lactones
    ]
    
    feature_match = False
    for pattern in feature_patterns:
        if mol.HasSubstructMatch(pattern):
            feature_match = True
            break
    
    if not feature_match:
        return False, "No common phenylpropanoid features found"
    
    return True, "Contains phenylpropane skeleton and common phenylpropanoid features"