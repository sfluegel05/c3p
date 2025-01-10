"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4_(smiles: str):
    """
    Classifies whether a chemical is a 3-hydroxy fatty acyl-CoA(4-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to match a hydroxyl group at the 3rd position of a fatty acyl chain
    hydroxyl_pattern = Chem.MolFromSmarts('CC(O)C(=O)')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3-hydroxy fatty acyl chain found"

    # Pattern to match CoA moiety, focusing on the presence of phosphate and nucleotide structure
    coa_pattern = Chem.MolFromSmarts('OP([O-])([O-])=O')
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if len(coa_matches) < 2:
        return False, "CoA moiety with proper phosphate groups not found"

    # Check for length of the carbon chain to be considered as a long chain typical of fatty acyl
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:  # Assuming a minimal chain length
        return False, "Carbon chain too short for fatty acyl"

    return True, "Matches 3-hydroxy fatty acyl-CoA(4-) structure"