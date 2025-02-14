"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:17683 1-acyl-sn-glycero-3-phosphoserine

An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
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
    
    # Look for phosphate group attached to the middle carbon of the glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("[CH2][CH1X4][CH2][OX2][PX4]([OX2])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found on the glycerol backbone"
    
    # Look for serine head group (more general pattern)
    serine_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[OX2H]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine head group found"
    
    # Look for acyl chain attached to the first carbon of the glycerol backbone
    acyl_pattern = Chem.MolFromSmarts("[CH2][OX2]C(=[OX1])[CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, "Incorrect number of acyl chains found"
    
    # Check acyl chain length (at least 6 carbons)
    acyl_chain_atoms = [mol.GetAtomWithIdx(idx) for idx in acyl_matches[0]]
    acyl_chain = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(*pair).GetIdx() for pair in zip(acyl_chain_atoms, acyl_chain_atoms[1:])])
    acyl_chain_length = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)
    if acyl_chain_length < 6:
        return False, "Acyl chain too short"
    
    return True, "Contains glycerol backbone with phosphate, serine head group, and acyl chain"