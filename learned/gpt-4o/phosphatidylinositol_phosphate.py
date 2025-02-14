"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

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

    # Check for glycerol backbone pattern: typically as O-C-C(O)-C(O)
    glycerol_pattern = Chem.MolFromSmarts("OC(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for inositol ring pattern: cyclohexane-1,2,3,4,5,6-hexol
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1(O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring structure found"

    # Check for the presence of one or more phosphate groups: -P(=O)(O)O-
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # Check for two long carbon chain esters: long acyl chains
    long_chain_pattern = Chem.MolFromSmarts("C(=O)OCCCCC")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Not enough long carbon chains"

    # Check that the molecular weight is within typical range but broaden criteria
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600 or mol_wt > 1500:
        return False, f"Molecular weight out of typical range for phosphatidylinositol phosphate: {mol_wt}"

    return True, "Contains glycerol backbone, inositol ring, requisite phosphate groups, and long ester-linked chains typical of phosphatidylinositol phosphates"