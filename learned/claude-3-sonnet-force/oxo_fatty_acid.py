"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:51748 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is a fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

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

    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Look for aldehydic (C=O) or ketonic (C(=O)C) groups
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C,H]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No aldehydic or ketonic groups found"

    # Check for carbon chain (at least 4 carbons)
    carbon_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    carbon_matches = mol.GetSubstructMatches(carbon_pattern)
    if not carbon_matches:
        return False, "Carbon chain too short (fewer than 4 carbons)"

    # Check molecular weight - fatty acids typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for fatty acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, "Too few carbons for fatty acid"
    if o_count < 2:
        return False, "Not enough oxygens (needs at least 2)"

    return True, "Contains carboxylic acid group and aldehydic/ketonic group within a fatty acid chain"