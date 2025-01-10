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
        str: Reason for the classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for inositol ring (6-membered ring with at least four hydroxy groups)
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"
    
    # Check for phosphorylated hydroxy groups on the inositol ring
    phosphorylated_inositol_pattern = Chem.MolFromSmarts("C(OP(=O)(O)O)C(O)C(O)C(O)C(O)C1")
    phosphorylated_matches = mol.GetSubstructMatches(phosphorylated_inositol_pattern)
    if not phosphorylated_matches:
        return False, "No phosphorylated inositol positions found"
    
    # Check for glycerol backbone linkage (OC[C@H] linked structures)
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone linkage found"
    
    # Check for the presence of two fatty acid chains (linked esters)
    fatty_acid_linkage_pattern = Chem.MolFromSmarts("OC(=O)C")
    links = len(mol.GetSubstructMatches(fatty_acid_linkage_pattern))
    if links < 2:
        return False, f"Found {links} ester linkages, require two for phosphoinositides"
        
    return True, "Contains inositol ring with specific phosphorylation and lipid linkages characteristic of a phosphoinositide"