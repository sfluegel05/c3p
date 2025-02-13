"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: CHEBI:18294 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerolipid with a phosphate group ester-linked to a carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]([OX2])[OX2][OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(~O)(~O)(~O)(~O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for ester linkage between phosphate and any carbon of the glycerol
    ester_pattern = Chem.MolFromSmarts("[OX2][CH2X4][CH2X4][CH2X4][OX2]P(~O)(~O)(~O)(~O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Phosphate not ester-linked to glycerol backbone"
    
    # Look for fatty acid chains (long carbon chains attached to oxygens of the glycerol)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[OX2][CH2X4][CH2X4][CH2X4]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chains attached to glycerol backbone"

    return True, "Contains glycerol backbone with phosphate group ester-linked and fatty acid chains"