"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Loose pattern for inositol ring (6-membered ring with at least 4 hydroxy groups)
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"
    
    # Check for any phosphate groups linked to oxygen (phosphodiester linkage)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphodiester found on inositol"

    # Check for glycerol-like linkage for lipid (O-C(=O)C side chain)
    lipid_linkage_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(lipid_linkage_pattern):
        return False, "No lipid linkage found"

    return True, "Contains inositol ring with phosphorylation and lipid linkages characteristic of a phosphoinositide"