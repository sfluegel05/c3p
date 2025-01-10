"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:35748 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (-CONH- or -CONR-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Check for a long hydrocarbon chain (more than 6 carbons)
    # We look for a chain of at least 6 carbons connected to the amide group
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons to be a fatty amide"

    # Check molecular weight - fatty amides typically have a higher molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a fatty amide"

    return True, "Contains an amide group with a long hydrocarbon chain"