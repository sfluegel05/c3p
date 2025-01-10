"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:35741 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three substituent groups 
    (acyl, alkyl, or alk-1-enyl) at positions sn-1, sn-2, and sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    
    # Look for ether groups (-O-C-)
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ether_matches = len(mol.GetSubstructMatches(ether_pattern)) - ester_matches  # Subtract esters as they also match ether pattern
    
    # Look for enol ether groups (-O-C=C-)
    enol_ether_pattern = Chem.MolFromSmarts("[OX2][CX3]=[CX3]")
    enol_ether_matches = len(mol.GetSubstructMatches(enol_ether_pattern))
    
    total_linkages = ester_matches + ether_matches + enol_ether_matches
    
    # Must have exactly 3 substituents total
    if total_linkages < 3:
        return False, f"Found only {total_linkages} substituent groups, need 3"
    
    # Check for long carbon chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = len(mol.GetSubstructMatches(chain_pattern))
    if chain_matches < 3:
        return False, "Missing long carbon chains"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short"

    # Check molecular weight - should be >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for triradylglycerol"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbons for triradylglycerol"

    # Count oxygens (minimum 3 for glycerol backbone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Must have at least 3 oxygens"

    # Success if we reach here
    substituent_types = []
    if ester_matches > 0:
        substituent_types.append("acyl")
    if ether_matches > 0:
        substituent_types.append("alkyl")
    if enol_ether_matches > 0:
        substituent_types.append("alk-1-enyl")
        
    return True, f"Contains glycerol backbone with {', '.join(substituent_types)} substituents"