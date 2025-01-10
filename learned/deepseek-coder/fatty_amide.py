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
    # More inclusive pattern that matches both primary and secondary amides
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,H1]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Check for a long hydrocarbon chain (at least 6 carbons total)
    # Count all carbons except those in the amide group
    c_count = sum(1 for atom in mol.GetAtoms() 
                 if atom.GetAtomicNum() == 6 and 
                 not atom.GetIdx() in [a.GetIdx() for a in mol.GetSubstructMatch(amide_pattern)])
    
    if c_count < 6:
        return False, f"Not enough carbons ({c_count}) for a fatty chain"

    # Check molecular weight - fatty amides typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a fatty amide"

    # Check for typical fatty acid characteristics:
    # At least 6 carbons in a chain (not necessarily consecutive)
    # Can include double/triple bonds and branches
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No suitable hydrocarbon chain found"

    # Additional check: at least 50% of carbons should be in chains
    # (excluding amide group carbons)
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count/total_c < 0.5:
        return False, "Too many carbons in non-chain positions"

    return True, "Contains an amide group with a long hydrocarbon chain"