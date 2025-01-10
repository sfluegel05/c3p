"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is a glycerophosphoinositol, reason for classification)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for inositol group - cyclohexane ring with 6 OH groups
    # More flexible pattern without specific stereochemistry
    inositol_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol group found"

    # Check for glycerol backbone with specific sn-3 phosphate attachment
    glycerol_phosphate = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]OP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "No glycerol backbone with phosphate at sn-3 position"

    # Check for phosphate-inositol connection
    # More specific pattern for phosphate connecting to inositol ring
    phosphoinositol_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC1C(O)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(phosphoinositol_pattern):
        return False, "Phosphate not properly connected to inositol"

    # Check for at least one ester group (fatty acid attachment)
    # More specific pattern for fatty acid esters
    ester_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]OP(=O)(O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester groups found in correct position"

    # Verify the complete core structure
    core_pattern = Chem.MolFromSmarts("[CH2X4]([CX4])[CHX4]([CX4])[CH2X4]OP(=O)(O)OC1C(O)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core glycerophosphoinositol structure not found"

    # Count key atoms and verify basic composition
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if o_count < 11:
        return False, f"Too few oxygen atoms ({o_count})"
    if c_count < 9:
        return False, f"Too few carbon atoms ({c_count})"
    if p_count != 1:
        return False, f"Must have exactly one phosphorus atom"

    # Additional check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]OC(=O)[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No proper fatty acid chain attachment found"

    return True, "Valid glycerophosphoinositol structure with correct core and substitution pattern"