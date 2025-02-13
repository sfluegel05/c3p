"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    Long-chain fatty acyl-CoA(4-) are characterized by a long fatty acid chain attached to a Coenzyme A structure
    with a net charge of 4- due to deprotonated phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)CO[P@]([O-])(=O)O[P@]([O-])(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A structure found"

    # Define pattern for an aliphatic tail with at least 12 carbons (flexible definition)
    # This assumes single bonds, but allows for variations typical in SMILES segments
    fatty_acyl_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")  # 12 linear CH2 groups
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No long-chain fatty acyl group detected"

    # Ensure deprotonated phosphates (di-phosphate pattern specific to these molecules)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-])([O-])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 4:
        return False, "Insufficient number of deprotonated phosphate groups"

    return True, "Classified as a long-chain fatty acyl-CoA(4-)"