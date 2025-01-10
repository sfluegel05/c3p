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

    # Update inositol detection with flexibility for hydroxyls
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol-like structure with hydroxyl groups found"

    # Adjusting phosphodiester linkage pattern to accommodate possible stereochemical degrees
    phosphodiester_pattern = Chem.MolFromSmarts("O[P](=O)([O-])O[C@H]")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage (O-P(=O)-O-C) found"

    # Check for ester linkages from glycerol to fatty acids
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond indicating fatty acid chains linked to glycerol found"

    # Check for long hydrocarbon chains using a more dynamic pattern by counting carbons
    carbon_chain_pattern = Chem.MolFromSmarts("[0*]-[0*]-[0*]-[0*]-[0*]-[0*]-[0*]-[0*]")
    c_chain_match = False
    for match in mol.GetSubstructMatches(carbon_chain_pattern):
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count >= 18:
            c_chain_match = True
            break

    if not c_chain_match:
        return False, "No sufficiently long aliphatic chains typical of fatty acids found"

    return True, "The molecule matches the structure of phosphatidylinositol: inositol ring with hydroxyls, phosphate linkage, glycerol esterified to fatty acid chains."