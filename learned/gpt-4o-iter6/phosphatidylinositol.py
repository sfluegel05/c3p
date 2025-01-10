"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has a phosphatidyl group esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Inositol recognition (six-membered ring with multiple hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("C1(CO)COC(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol-like structure with hydroxyl groups found"

    # Detect phosphodiester linkage (glycerophosphate connected to inositol).
    # Consider variations of -OPO ester bond linking inositol and glycerol
    phosphodiester_pattern = Chem.MolFromSmarts("O[P](=O)([O-])O[C@@H](C)O")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage (O-P(=O)-O-C) found to inositol"

    # Check for at least one ester linkage typical of fatty acids attached to glycerol
    fatty_acid_ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CO)CO")
    ester_matches = mol.GetSubstructMatches(fatty_acid_ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester bond indicating fatty acid chains linked to glycerol found"

    # Check for long hydrocarbon chains typical of phosphatidyl groups (18+ carbons)
    long_chain_pattern = Chem.MolFromSmarts("C-C-C-C-C-C-C-C")  # Adjust to match long chains
    if not any(len(Chem.rdmolops.GetMolFrags(Chem.DeleteSubstructs(mol, match, True, True), False, False)) >= 18 for match in mol.GetSubstructMatches(long_chain_pattern)):
        return False, "No sufficiently long aliphatic chains typical of fatty acids found"

    return True, "Matches the structure of phosphatidylinositol: inositol linked via a phosphate ester to glycerol containing long chain fatty acids"