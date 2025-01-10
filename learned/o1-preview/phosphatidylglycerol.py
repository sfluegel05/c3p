"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    Phosphatidylglycerols consist of a glycerol backbone with fatty acyl chains at sn-1 and sn-2 positions,
    and a phosphatidyl group attached to another glycerol at the sn-3 position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove stereochemistry to generalize matching
    Chem.RemoveStereochemistry(mol)
    
    # Pattern for glycerol backbone with two ester-linked fatty acyl chains
    esterified_glycerol = Chem.MolFromSmarts("COC(=O)[#6]")
    ester_matches = mol.GetSubstructMatches(esterified_glycerol)
    if len(ester_matches) < 2:
        return False, f"Less than two ester-linked fatty acyl chains found ({len(ester_matches)} found)"

    # Pattern for phosphate group attached to glycerol (phosphatidyl group)
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(OCCO)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphatidylglycerol headgroup found"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(CO)(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Passed all checks
    return True, "Molecule matches phosphatidylglycerol structure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'phosphatidylglycerol',
        'definition': 'A glycerophosphoglycerol that is glycerol in which the hydrogen of one of the primary hydroxy groups has been replaced by a phosphatidyl group.',
    },
}