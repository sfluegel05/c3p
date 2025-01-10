"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and typically an aldehyde functional group or exists as a cyclic form.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 5 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"

    # Check for aldehyde group: pattern C(=O) connected to a terminal H (open-chain form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH1]")
    # Also check for cyclic ring forms (furanose/pyranose)
    ring_info = mol.GetRingInfo()
    has_furanose_pyranose = any(len(ring) in [5, 6] for ring in ring_info.AtomRings())
    
    # Check for undesired groups like phosphates or carboxylic acids
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, not a simple aldopentose"
    if mol.HasSubstructMatch(carboxy_pattern):
        return False, "Contains carboxyl group, not a simple aldopentose"

    if not mol.HasSubstructMatch(aldehyde_pattern) and not has_furanose_pyranose:
        return False, "No aldehyde group or furanose/pyranose ring detected"

    # Check for sufficient oxygen atoms (4 minimum for the hydroxyl groups on the carbons)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Insufficient oxygen atoms: expected at least 4, found {o_count}"

    return True, "Molecule matches the characteristics of an aldopentose"