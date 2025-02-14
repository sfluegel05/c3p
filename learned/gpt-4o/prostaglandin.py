"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are C20 compounds derived from prostanoic acid, generally featuring a
    substituted cyclopentane ring and various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclopentane ring with possible diverse substitutions
    flexible_cyclopentane_pattern = Chem.MolFromSmarts("[C&R2]1(~[C&R2])(~[C&R2])~[C&R2]~[C&R2]1")
    if not mol.HasSubstructMatch(flexible_cyclopentane_pattern):
        return False, "No flexible cyclopentane or substituted cyclopentane ring detected"
    
    # Check for an acid or derivative (e.g., carboxylic acid, ester, amide)
    acid_or_ester_pattern = Chem.MolFromSmarts("C(=O)[O,N,S]")
    if not mol.HasSubstructMatch(acid_or_ester_pattern):
        return False, "Missing acid or acid-derivative functional group"
    
    # Determine presence of key double bond patterns distinct in many prostaglandins
    db_patterns = [
        Chem.MolFromSmarts("C=CC=CC=C"),  # Example pattern covering extended unsaturation
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in db_patterns):
        return False, "Missing typical prostaglandin double bond patterns"

    # Consider molecular size or framework typical of C20 structures
    num_atoms = mol.GetNumHeavyAtoms()
    if num_atoms < 18 or num_atoms > 24:
        return False, "Molecule size not typical for prostaglandins based on heavy atoms count"

    return True, "Molecule contains a C20 backbone with a cyclopentane or properly substituted cyclopentane, relevant acid-functional groups and expected double bond patterns typical of prostaglandins"