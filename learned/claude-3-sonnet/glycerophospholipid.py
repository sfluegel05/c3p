"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid has a glycerol backbone with a phosphate group ester-linked 
    to a terminal carbon and at least one fatty acid chain.

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

    # Check for phosphate group (-P(=O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([O,N])[O,N]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CH2X4][CHX4][CH2X4][OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Verify phosphate is connected to glycerol
    phospho_glycerol = Chem.MolFromSmarts("[OX2][CH2X4][CHX4][CH2X4][OX2]P(=O)([O,N])[O,N]")
    if not mol.HasSubstructMatch(phospho_glycerol):
        return False, "Phosphate not connected to glycerol backbone"

    # Look for at least one ester group (-O-C(=O)-)
    # Note: Some glycerophospholipids can have ether linkages instead
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    
    if not (has_ester or has_ether):
        return False, "No ester or ether linkages found"

    # If we find esters, check for fatty acid chains
    if has_ester:
        fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        if not mol.HasSubstructMatch(fatty_acid_pattern):
            return False, "No fatty acid chains found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 8:
        return False, "Too few carbons for glycerophospholipid"
    if o_count < 4:
        return False, "Too few oxygens for glycerophospholipid"
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"

    # Verify molecular weight - glycerophospholipids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for glycerophospholipid"

    return True, "Contains glycerol backbone with phosphate group and appropriate fatty acid/ether chains"