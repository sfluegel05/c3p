"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid has a glycerol backbone, one acyl chain, and a phosphate group,
    with no additional phospholipid head groups like choline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone - it should include an ether group linking to a central carbon
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for the phosphate group; phosphoric acid ester
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for monoacyl group (-O-C(=O)-R) - only one such attachment should be present
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"

    # Ensure no additional head groups like choline are present
    choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C(O)1")
    if mol.HasSubstructMatch(choline_pattern) or mol.HasSubstructMatch(inositol_pattern):
        return False, "Found head group(s) such as choline or inositol, not a lysophosphatidic acid"

    # Additional checkpoint: Verify there are no extra phospholipid features
    num_phosphoryl_groups = Descriptors.CalcNumRotatableBonds(mol)
    if num_phosphoryl_groups > 1:
        return False, "Complex phospholipid structure detected, excessive rotatable bonds"

    return True, "Contains glycerol backbone, one acyl group, and a phosphate group without additional head groups"