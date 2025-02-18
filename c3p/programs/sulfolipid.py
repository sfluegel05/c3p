"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:73404 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfate group (O-SO3H) attached to a lipid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfate group (O-SO3H or O-SO3-)
    # SMARTS pattern for sulfate ester (O-S(=O)(=O)[O-] or O-S(=O)(=O)O)
    sulfate_pattern = Chem.MolFromSmarts("[O][S](=O)(=O)[O]")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group (O-SO3H) found"

    # Check lipid characteristics: significant carbon count and molecular weight
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Insufficient carbons ({c_count}) for lipid structure"

    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for lipid"

    # Additional check for long-chain fatty acid ester/amide
    # Look for at least one ester or amide group connected to a chain of >=12 carbons
    # Using a simple pattern for long aliphatic chains (>=12 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long aliphatic chain (>=12 carbons) detected"

    return True, "Contains sulfate group attached to lipid structure with long chains"