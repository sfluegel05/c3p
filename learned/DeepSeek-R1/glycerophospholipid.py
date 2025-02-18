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

    # Check for glycerol backbone (three carbons connected as C-O-C-O-C)
    glycerol = Chem.MolFromSmarts("[CH2]-[CH](-[OX2])-[CH2]")
    if not mol.HasSubstructMatch(glycerol):
        return False, "No glycerol backbone"

    # Find phosphate group (P connected to four oxygens, at least one ester linkage)
    phosphate = Chem.MolFromSmarts("[PX4](=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate)
    if not phosphate_matches:
        return False, "No phosphate group"

    # Verify phosphate is ester-linked to glycerol's terminal carbon
    # Check if the phosphate is connected via an oxygen to the central carbon of glycerol
    ester_link = Chem.MolFromSmarts("[CH2]-[CH](-[OX2]-P(=O)([OX2])[OX2])-[CH2]")
    if not mol.HasSubstructMatch(ester_link):
        return False, "Phosphate not ester-linked to glycerol"

    # Check for two ester groups (fatty acids) on the other two oxygens of glycerol
    ester = Chem.MolFromSmarts("[CX4][OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester)
    if len(ester_matches) < 2:
        return False, f"Only {len(ester_matches)} ester groups, need at least 2"

    return True, "Glycerol backbone with phosphate and two ester groups"