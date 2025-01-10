"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and an aldehyde functional group.

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
    
    # Check for aldehyde group: pattern C=O connected to an H and a carbon
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH1]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        # If no aldehyde is found, it might be in cyclic form, check for the ring
        ring_info = mol.GetRingInfo()
        has_furanose_pyranose = any(ring_size in [5, 6] for ring_size in ring_info.AtomRings())
        if not has_furanose_pyranose:
            return False, "No aldehyde group or furanose/pyranose ring detected"
    
    # Check for oxygen atoms - each carbon in aldopentoses should be attached to an oxygen
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Insufficient oxygen atoms: expected at least 4, found {o_count}"
    
    return True, "Molecule matches the characteristics of an aldopentose"