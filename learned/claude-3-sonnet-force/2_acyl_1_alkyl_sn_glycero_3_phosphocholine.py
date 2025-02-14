"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:18035 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule belongs to the class 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
    based on its SMILES string. This class is defined as an alkyl,acyl-sn-glycero-3-phosphocholine
    with unspecified alkyl and acyl groups located at positions 1 and 2 respectively.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
    
    # Get the glycerol atoms
    glycerol_atoms = mol.GetSubstructMatch(glycerol_pattern)
    
    # Check for acyl group at position 2
    acyl_group = Chem.MolFromSmiles("C(=O)") # Simple acyl group
    acyl_match = mol.GetSubstructMatches(acyl_group, maxMatches=1, uniqueOnly=True)
    if not acyl_match:
        return False, "No acyl group found at position 2"
    acyl_idx = acyl_match[0][0]
    if acyl_idx != glycerol_atoms[1]:
        return False, "Acyl group not at position 2"
    
    # Check for alkyl group at position 1
    alkyl_atoms = [atom.GetIdx() for atom in mol.GetAtomWithIdx(glycerol_atoms[0]).GetNeighbors() if atom.GetAtomicNum() == 6]
    if not alkyl_atoms:
        return False, "No alkyl group found at position 1"
    
    # Count carbon atoms in the alkyl chain
    alkyl_chain = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCC") # Assume max 18 carbon atoms
    alkyl_match = mol.GetSubstructMatches(alkyl_chain, maxMatches=1, uniqueOnly=True)
    if alkyl_match:
        alkyl_start_idx = alkyl_match[0][0]
        alkyl_length = len([atom.GetIdx() for atom in mol.GetAtomWithIdx(alkyl_start_idx).GetNeighbors() if atom.GetAtomicNum() == 6])
        if alkyl_length < 4:
            return False, "Alkyl chain too short"
    else:
        return False, "Alkyl chain not found"
    
    return True, "Molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"