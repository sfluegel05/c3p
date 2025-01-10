"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: CHEBI:73009 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    A 1-O-acylglycerophosphoethanolamine has a glycerol backbone with a phosphoethanolamine group at the 3-position
    and an acyl group at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with specific connectivity
    glycerol_pattern = Chem.MolFromSmarts("[*:1][C@H]([OH])[C@H]([OH])[*:2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No proper glycerol backbone found"

    # Look for phosphoethanolamine group at the 3-position (including protonated forms)
    phosphoethanolamine_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])OCC[NH2,NH3+]")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found at the 3-position"

    # Look for acyl group specifically at the 1-position
    # Pattern: [C@H](O)COC(=O) for sn-1 position
    acyl_pattern = Chem.MolFromSmarts("[C@H]([OH])COC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group found at the 1-position"

    # Verify that the acyl group is connected to a long chain
    # Look for at least 10 carbons in the acyl chain
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Acyl chain too short"

    # Count total carbons to ensure sufficient size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbons for 1-O-acylglycerophosphoethanolamine"

    return True, "Contains glycerol backbone with phosphoethanolamine at 3-position and acyl group at 1-position"