"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate typically contains a glycerol backbone, an inositol ring, multiple phosphate groups, 
    and ester-linked long fatty acid chains.

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

    # Check for glycerol backbone pattern: typically as -C(O)C(X)CO- where X is variable
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for inositol ring pattern, generalize pattern for hydroxylated cyclohexane
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1(O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring structure found"

    # Check for at least one phosphate group: -P(=O)(O)O-
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # More adaptable pattern for long carbon chains, default to a length >=8
    long_chain_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCC")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Not enough long carbon chains, possibly due to incorrect pattern"

    # Check molecular weight to confirm large molecule
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600 or mol_wt > 1500:
        return False, f"Molecular weight out of typical range for phosphatidylinositol phosphate: {mol_wt}"

    return True, "Contains glycerol backbone, inositol ring, requisite phosphate groups, and long ester-linked chains typical of phosphatidylinositol phosphates"