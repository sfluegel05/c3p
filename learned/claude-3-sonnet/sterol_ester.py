"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More flexible steroid core pattern that matches different steroid variants
    # Matches the basic 6-6-6-5 ring system with flexibility in bond types
    steroid_core = Chem.MolFromSmarts("C1C[C,=C]2[C,=C][C,=C][C,=C]3[C,=C][C,=C][C,=C]4[C,=C][C,=C][C,=C]4[C,=C]3[C,=C][C,=C]2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # More flexible pattern for sterol ester at position 3
    # This pattern allows for different steroid core variants while maintaining the key ester position
    sterol_ester_pattern = Chem.MolFromSmarts("[C,=C]1[CH2][CH2,CH][C,=C]([OX2][CX3](=[OX1])[#6])[C,=C]2[C,=C][C,=C][C,=C]3[C,=C][C,=C][C,=C]4[C,=C][C,=C][C,=C]4[C,=C]3[C,=C]2[C,=C]1")
    if not mol.HasSubstructMatch(sterol_ester_pattern):
        return False, "Ester group not at correct position on steroid core"

    # Count carbons (sterols typically have >20 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for a sterol ester"
    
    # Count rings (sterols typically have 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Too few rings for a sterol structure"
    
    # Look for characteristic branching that's common in sterols
    # This helps distinguish sterols from other steroid types
    methyl_branch = Chem.MolFromSmarts("[CH3][C]([C,=C])([C,=C])[C,=C]")
    if not mol.HasSubstructMatch(methyl_branch):
        return False, "Missing characteristic sterol branching"
        
    return True, "Contains steroid core with ester group at correct position"