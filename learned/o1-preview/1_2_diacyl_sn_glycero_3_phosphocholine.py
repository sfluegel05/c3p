"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 1,2-diacyl-sn-glycero-3-phosphocholine
Definition: The conjugate base of a 1,2-diacyl-sn-glycero-3-phosphocholine compound formed by deprotonation of the phosphate OH group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.
    This class of compounds consists of a glycerol backbone with fatty acid chains esterified at the
    sn-1 and sn-2 positions, and a phosphocholine group attached at the sn-3 position via a phosphate ester linkage.
    The molecule is in its deprotonated form at the phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has been sanitized
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.KekulizeException:
        return False, "Unable to sanitize molecule"

    # Define SMARTS patterns

    # Pattern for glycerol backbone with two esterified positions and one phosphate ester
    glycerol_pattern = Chem.MolFromSmarts("C(O[CX3](=O)[C,c])(O[CX3](=O)[C,c])CO")

    # Match the glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with two ester groups not found"

    # Pattern for phosphocholine group with deprotonated phosphate
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine group not found"

    # Check for two esterified fatty acid chains at sn-1 and sn-2 positions
    ester_pattern = Chem.MolFromSmarts("OC(=O)[C,c]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Less than two esterified fatty acid chains found"

    # Check for deprotonated phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group is not deprotonated"

    return True, "Molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine with deprotonated phosphate group"

__metadata__ = {   
    'chemical_class': {   
        'name': '1,2-diacyl-sn-glycero-3-phosphocholine',
        'definition': 'The conjugate base of a 1,2-diacyl-sn-glycero-3-phosphocholine compound formed by deprotonation of the phosphate OH group.',
    },
}