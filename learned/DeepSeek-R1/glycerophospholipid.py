"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: CHEBI:28874 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid has a glycerol backbone with a phosphate group ester-linked to a terminal carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for glycerol backbone (three carbons in a row with oxygen substituents)
    # Pattern allows for any substituents (including esters) on the carbons
    glycerol = Chem.MolFromSmarts("[C;H2,H1]([OX2])-[C;H1]([OX2])-[C;H2,H1]([OX2])")
    if not mol.HasSubstructMatch(glycerol):
        return False, "No glycerol backbone"

    # Find phosphate group connected via ester to a terminal carbon of glycerol
    # Pattern: terminal CH2-O-P (phosphate ester)
    phosphate_ester = Chem.MolFromSmarts("[CH2]-[OX2]-P(=O)([OX2])-[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_ester)
    if not phosphate_matches:
        return False, "No phosphate ester group on terminal carbon"

    # Verify the phosphate is attached to glycerol's terminal carbon
    # Check if any phosphate ester oxygen is connected to glycerol's terminal carbon
    glycerol_phosphate = Chem.MolFromSmarts(
        "[C;H2,H1]([OX2])-[C;H1]([OX2])-[C;H2,H1]([OX2]-P(=O)([OX2])-[OX2])"
    )
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "Phosphate not ester-linked to glycerol's terminal carbon"

    # Check for at least one fatty acid ester (long chain)
    ester_pattern = Chem.MolFromSmarts("[CX4][OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No fatty acid ester groups detected"

    # Optional: Check for typical chain length in fatty acids (e.g., at least 8 carbons)
    # This is heuristic and may vary, but helps exclude short-chain esters
    chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "Insufficient chain length for fatty acids"

    return True, "Glycerol backbone with terminal phosphate ester and fatty acid(s)"