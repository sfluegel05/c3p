"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol with one or more phosphate groups on the inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        (bool, str): Tuple where first element indicates if the molecule is a phosphoinositide,
                      and the second element provides the reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for potential glycerol backbone with fatty acids (ester connections for flexibility in patterns)
    glycerol_patterns = [
        Chem.MolFromSmarts("[CX3](=O)OC[CH](O)"),
        Chem.MolFromSmarts("OCC(OC(=O))")  # Allowing different orientations
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns):
        return False, "No glycerol backbone with fatty acids found"

    # Inositol ring detection: Stereochemistry may vary or not be depicted
    inositol_patterns = [
        Chem.MolFromSmarts("C1([C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)O)"),
        Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C1O)")  # Non-chiral
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in inositol_patterns):
        return False, "Inositol ring not found"

    # Check for phosphate groups bonded to the inositol ring
    phosphate_patterns = [
        Chem.MolFromSmarts("C(OP(=O)(O)O)C"),
        Chem.MolFromSmarts("C(OP(O)(O)=O)")  # Different pattern recognition
    ]

    phosphate_matches = [mol.GetSubstructMatches(pattern) for pattern in phosphate_patterns]

    # Ensure at least one phosphate attachment is found
    if not any(len(matches) > 0 for matches in phosphate_matches):
        return False, "No phosphorylated inositol ring found"

    return True, "Molecule has a phosphatidylinositol backbone with one or more phosphate groups on the inositol ring"