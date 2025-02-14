"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate contains a glycerol backbone, an inositol ring,
    multiple phosphate groups, and ester-linked long fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone pattern installed with flexibility
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)C(O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for inositol ring pattern with flexible stereochemistry
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring structure found"

    # Look for phosphate groups associated with inositol
    phosphate_matches = Chem.MolFromSmarts("OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_matches):
        return False, "No phosphate groups found"

    # Recognize structures with typical ester-linked fatty acid chains (possibly 14-20 carbons)
    ester_chain_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(ester_chain_pattern):
        return False, "Not enough ester-linked long carbon chains"

    # Validate using molecular weight if other criteria validate but molecular identity ambiguities persist
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600 or mol_wt > 1600:
        return False, f"Molecular weight out of typical range for phosphatidylinositol phosphate: {mol_wt}"

    return True, "Contains attributes consistent with phosphatidylinositol phosphate: glycerol backbone, inositol ring, requisite phosphate groups, and long ester-linked fatty acid chains"