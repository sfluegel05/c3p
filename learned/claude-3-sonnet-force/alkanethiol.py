"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:25710 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for -SH group
    sh_pattern = Chem.MolFromSmarts("[SH]")
    if not mol.HasSubstructMatch(sh_pattern):
        return False, "No sulfanyl (-SH) group found"
    
    # Look for linear alkyl chains attached to sulfur (avoiding ring structures)
    alkyl_pattern = Chem.MolFromSmarts("[SH]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No linear alkyl chain attached to sulfur"

    # Allow additional functional groups like -OH, -C=O, -P=O, etc.
    # but exclude sulfur-containing ring structures or additional substituted sulfur atoms
    valid_pattern = Chem.MolFromSmarts("[SH]~[CX4]~[CX4]~[!#16,#6,#1]")
    if not mol.HasSubstructMatch(valid_pattern):
        return False, "Contains additional sulfur-containing groups or ring structures"

    return True, "Contains sulfanyl (-SH) group attached to a linear alkyl chain"