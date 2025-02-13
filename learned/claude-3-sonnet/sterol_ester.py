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

    # Basic steroid core pattern - four fused rings (three 6-membered, one 5-membered)
    # More general pattern that doesn't specify aromaticity
    steroid_core = Chem.MolFromSmarts("C1C2CCC3CCCC4CCCC(C1)C24C3")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatches(ester_pattern):
        return False, "No ester group found"

    # Pattern for 3-position oxygen (more general pattern)
    # This matches the oxygen at position 3 of the first ring
    sterol_3_o = Chem.MolFromSmarts("C1CC([OX2])CC2C1")
    if not mol.HasSubstructMatch(sterol_3_o):
        return False, "No oxygen at position 3 of steroid core"

    # Count carbons (allowing for smaller sterol esters)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Lowered threshold
        return False, f"Too few carbons ({c_count}) for a sterol ester"

    # Count rings (must have at least 4 for steroid core)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Too few rings for a sterol structure"

    # Look for at least one angular methyl group (more general pattern)
    angular_methyl = Chem.MolFromSmarts("[CH3][C]([#6])([#6])[#6]")
    if not mol.HasSubstructMatch(angular_methyl):
        return False, "Missing characteristic angular methyl group"

    # Additional check for connectivity between ester and steroid core
    # The oxygen of the ester must be connected to the steroid core
    ester_o = Chem.MolFromSmarts("[OX2]([#6])[CX3](=[OX1])[#6]")
    matches = mol.GetSubstructMatches(ester_o)
    if not any(mol.GetRingInfo().IsAtomInRingOfSize(match[0], 6) for match in matches):
        return False, "Ester group not properly connected to steroid core"

    return True, "Contains steroid core with ester at position 3 and characteristic features"