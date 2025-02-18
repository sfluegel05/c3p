"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: CHEBI:17400 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is formed by esterification of cholesterol's 3-hydroxy group with a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define cholesterol core pattern (approximate structure)
    # Looking for the sterol skeleton with ester oxygen at position 3
    cholesterol_core = Chem.MolFromSmarts("[C@H]1([C@H](C[C@@H]2[C@@]1(C3=C(C[C@H](OC(=O)*)C3)CC2)C)[H])CC[C@@H]4[C@@]5(CC[C@H](C4)CC5)C)")
    
    if not mol.HasSubstructMatch(cholesterol_core):
        return False, "Cholesterol core not found or ester group missing at position 3"

    # Check for exactly one ester group connected to the core
    ester_pattern = Chem.MolFromSmarts("[OX2][C]=O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester group, found {len(ester_matches)}"

    # Verify the ester is connected to a carbon chain (fatty acid)
    # Look for at least 2 carbons in the acid part (since even unsaturated acids have chains)
    acid_chain = Chem.MolFromSmarts("[CX4,CX3]=O.[OX2][C]=O")
    if not mol.HasSubstructMatch(acid_chain):
        return False, "Ester does not connect to a sufficient carbon chain"

    return True, "Cholesterol ester with 3-O-acyl group"