"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem

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
    
    # Updated pattern for cyclopentane ring allowing more variation
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCC(C1)O")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No flexible cyclopentane or substituted cyclopentane ring detected"
    
    # Check for an acid or derivative functional group (e.g., carboxylic acid, ester, amide)
    acid_or_ester_pattern = Chem.MolFromSmarts("COC=O | C(=O)[O,N,S]")
    if not mol.HasSubstructMatch(acid_or_ester_pattern):
        return False, "Missing typical acid or acid-derivative functional group"
    
    # Add alternate double bond pattern considerations more aligned with prostaglandin molecules
    db_pattern_subtle = Chem.MolFromSmarts("C/C=C/C")
    if not mol.HasSubstructMatch(db_pattern_subtle):
        return False, "Missing expected double bond patterns in the side chain"

    # Consider molecular size or framework typical of C20 structures
    num_atoms = mol.GetNumHeavyAtoms()
    if num_atoms < 18 or num_atoms > 24:
        return False, "Molecule size not typical for prostaglandins based on heavy atoms count"

    return True, "Molecule contains essential features of C20 prostaglandins, including cyclopentane, functional groups, and double bond patterns."