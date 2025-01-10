"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: CHEBI:36044 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroperoxy group (-OOH) attached to a carbon
    hydroperoxy_pattern = Chem.MolFromSmarts("[CX4][OX2][OX1]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # Check for lipid characteristics (long carbon chain and carboxylic acid/ester)
    # Look for at least 10 carbons in a chain
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"

    # Check for carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or ester group found"

    # Check molecular weight - lipids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for lipid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for lipid"
    if o_count < 3:
        return False, "Too few oxygens for lipid hydroperoxide"

    return True, "Contains hydroperoxy group and lipid characteristics"