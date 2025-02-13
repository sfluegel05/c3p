"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A molecule is classified as such if it contains a glycerol backbone with specific stereochemistry,
    acyl substituents at the sn-1 and sn-2 positions, and a phosphoserine at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycerol backbone with correct stereochemistry
    glycerol_sn3_pattern = Chem.MolFromSmarts("[C@H](C([O])(=O)O)[O]C(C(=O)[*])[C@@H](O)COP=O")
    if not mol.HasSubstructMatch(glycerol_sn3_pattern):
        return False, "No proper glycerol sn-3 stereochemistry found"

    # Look for acyl groups on sn-1 and sn-2 positions
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H]C(=O)[*]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 2:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 2"

    # Look for phosphoserine moiety
    phosphoserine_pattern = Chem.MolFromSmarts("COP(=O)(O)OC[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    return True, "Contains the structures indicative of a 3-sn-phosphatidyl-L-serine"

# Example of metadata, not necessary for the function itself
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17855',
        'name': '3-sn-phosphatidyl-L-serine',
        'definition': 'A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.',
        'parents': ['CHEBI:47778', 'CHEBI:76886']
    }
}