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

    # Check for inositol ring (cyclohexane with 6 OH groups)
    # More flexible pattern that matches different stereoisomers
    inositol_pattern = Chem.MolFromSmarts("[OX2][CH1]1[CH1]([OX2])[CH1]([OX2])[CH1]([OX2])[CH1]([OX2])[CH1]1[OX2]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphate group connected to inositol
    # Pattern matches phosphate that bridges inositol and glycerol
    phosphate_pattern = Chem.MolFromSmarts("[OX2][CH1]-1-[CH1]([OX2])[CH1]([OX2])[CH1]([OX2])[CH1]([OX2])[CH1]-1[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group properly connected to inositol"

    # Check for glycerol backbone with correct connectivity
    # Pattern matches the glycerol backbone connected to phosphate
    glycerol_pattern = Chem.MolFromSmarts("[OX2]CC([OX2]C)COP(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct phosphate attachment"

    # Check for at least one ester group (fatty acid attachment)
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)OC[CH1]([OX2])CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bonds found for fatty acid attachment"

    # Verify basic composition
    atoms = mol.GetAtoms()
    p_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    
    if p_count != 1:
        return False, f"Must have exactly one phosphorus atom, found {p_count}"
    if o_count < 11:
        return False, f"Insufficient oxygen atoms for glycerophosphoinositol structure (minimum 11 required, found {o_count})"
    if c_count < 9:
        return False, f"Insufficient carbon atoms for glycerophosphoinositol structure (minimum 9 required, found {c_count})"

    # Verify complete core structure with proper connectivity
    core_pattern = Chem.MolFromSmarts("[OX2]CC([OX2]C)COP(=O)([OX2])[OX2][CH1]1[CH1]([OX2])[CH1]([OX2])[CH1]([OX2])[CH1]([OX2])[CH1]1[OX2]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core glycerophosphoinositol structure not found with proper connectivity"

    # Check for long carbon chains (fatty acids)
    carbon_chain = Chem.MolFromSmarts("CCCCC")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chains found for fatty acids"

    return True, "Valid glycerophosphoinositol structure with correct core components and substitution pattern"