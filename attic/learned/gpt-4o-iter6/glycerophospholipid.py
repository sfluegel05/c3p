"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerolipid with a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a glycerol backbone pattern allowing esterified carbons [OCC(O)CO]
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for a phosphate group ester-linked [P(=O)(O)O]
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group ester-linked to glycerol backbone"

    # Ensure there are at least two ester linkages, typical of fatty acid chains [C(=O)OCC]
    ester_pattern = Chem.MolFromSmarts("C(=O)OCC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2 for fatty acid chains"

    # Check molecular weight within a reasonable range for glycerophospholipids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (500 <= mol_wt <= 1000):
        return False, "Molecular weight outside typical glycerophospholipid range"

    return True, "Valid glycerophospholipid structure with glycerol backbone, phosphate group, and fatty acid chains"