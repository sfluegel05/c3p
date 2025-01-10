"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26195 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is aromatic
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule is not aromatic"

    # Look for phenylpropane skeleton pattern (C6H5-C3H7)
    phenylpropane_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(phenylpropane_pattern):
        return False, "No phenylpropane skeleton found"

    # Check for additional functional groups or modifications
    # Phenylpropanoids can have various functional groups, so we don't enforce strict rules here
    # but we can check for common ones like hydroxyl, methoxy, etc.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    # Count the number of carbons and hydrogens to ensure it fits the general formula
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

    # Phenylpropanoids typically have at least 9 carbons (C6H5-C3H7)
    if c_count < 9:
        return False, "Too few carbons for phenylpropanoid"

    # Check molecular weight - phenylpropanoids typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for phenylpropanoid"

    # If all checks pass, classify as phenylpropanoid
    return True, "Contains phenylpropane skeleton and is aromatic"