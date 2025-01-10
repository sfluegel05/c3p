"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid containing a trimethyl and flexible 17-furanyl or its derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for general triterpenoid skeleton (around 30 carbons, allowance due to variations)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27 or c_count > 35:  # Provide allowance for derivatizations
        return False, f"Number of carbon atoms {c_count} not expected for a triterpenoid"

    # Check for high oxygen content (limonoids are highly oxygenated)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Low oxygen count for highly oxygenated limonoid: {o_count}"

    # Flexibly consider characteristic oxygen functional groups
    ketone_pattern = Chem.MolFromSmarts("C=O")
    lactone_pattern = Chem.MolFromSmarts("C(=O)O")
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl groups
    if not (mol.HasSubstructMatch(ketone_pattern) or 
            mol.HasSubstructMatch(lactone_pattern) or
            mol.HasSubstructMatch(alcohol_pattern)):
        return False, "Lack of characteristic oxygenated functionalities (ketone, lactone, alcohol)"
    
    # Optionally check for generic trimethyl or equivalent functional groups
    trimethyl_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "Lack of characteristic trimethyl or equivalent groups for limonoid structure"

    return True, "SMILES string corresponds to a limonoid structure"