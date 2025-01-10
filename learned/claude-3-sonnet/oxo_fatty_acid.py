"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:35667 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid contains at least one aldehydic or ketonic group in addition 
    to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches > 1:
        return False, "Multiple carboxylic acid groups found"

    # Look for ketone groups - exclude amides and esters
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])([#6])[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Look for aldehyde groups
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    total_oxo_groups = len(ketone_matches) + len(aldehyde_matches)
    
    if total_oxo_groups == 0:
        return False, "No ketone or aldehyde group found"

    # Count atoms and check ratios
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Typical fatty acids should have C4-C30
    if c_count < 4 or c_count > 30:
        return False, "Carbon count outside typical fatty acid range (4-30)"
    
    # Check for nitrogen-containing compounds
    if n_count > 0:
        return False, "Contains nitrogen atoms"
    
    # Check oxygen to carbon ratio - should be relatively low for fatty acids
    if o_count / c_count > 0.5:
        return False, "Too many oxygen atoms for a fatty acid"

    # Look for sugar patterns
    sugar_pattern = Chem.MolFromSmarts("[OX2H1][CX4H1]([OX2])[CX4H1]([OX2])")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moiety"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 1:  # Allow only 1 ring max
        return False, "Too many rings for a fatty acid"

    # Look for long carbon chain characteristic of fatty acids
    long_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No sufficiently long carbon chain found"

    # Count rotatable bonds to verify chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Too rigid to be a fatty acid"

    # Check for excessive branching
    branching_pattern = Chem.MolFromSmarts("[CX4]([#6])([#6])([#6])[#6]")
    if len(mol.GetSubstructMatches(branching_pattern)) > 2:
        return False, "Too much branching for a fatty acid"

    # Success case - found carboxylic acid and at least one oxo group
    oxo_description = []
    if len(ketone_matches) > 0:
        oxo_description.append(f"{len(ketone_matches)} ketone(s)")
    if len(aldehyde_matches) > 0:
        oxo_description.append(f"{len(aldehyde_matches)} aldehyde(s)")
    
    return True, f"Contains carboxylic acid and {' and '.join(oxo_description)}"