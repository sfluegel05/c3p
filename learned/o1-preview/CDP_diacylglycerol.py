"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is a glycerol backbone with two fatty acid chains attached via ester bonds
    at positions 1 and 2, and a cytidine diphosphate group attached via a phosphodiester bond
    at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns without stereochemistry

    # Glycerol backbone: three carbons with two hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Ester linkage pattern: ester group
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;!$(C=O)]")  # Ester linkage
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Phosphodiester linkage pattern: phosphate group connected to glycerol
    phospho_glycerol_pattern = Chem.MolFromSmarts("O[P](=O)(O)OCCO")  # Phosphodiester linked to glycerol
    if not mol.HasSubstructMatch(phospho_glycerol_pattern):
        return False, "No phosphodiester linkage to glycerol found"

    # Cytidine moiety pattern: cytosine base connected to ribose sugar
    cytidine_pattern = Chem.MolFromSmarts("n1cccnc1O[C@H]2[C@@H](O)[C@H](O)[C@@H](CO[P](=O)(O)O[P](=O)(O)O)O2")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine moiety found"

    # Check that the phosphate groups are connected
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate groups not properly connected"

    # Ensure connectivity between glycerol backbone and ester groups
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    ester_atoms = set()
    for match in ester_matches:
        # The oxygen atom in the ester linkage connected to glycerol
        ester_oxygen_idx = match[2]
        ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
        for neighbor in ester_oxygen.GetNeighbors():
            if neighbor.GetIdx() in sum(glycerol_matches, ()):
                ester_atoms.add(ester_oxygen_idx)
    if len(ester_atoms) < 2:
        return False, "Ester groups not properly attached to glycerol backbone"

    # Ensure connectivity between glycerol backbone and phosphodiester group
    phospho_matches = mol.GetSubstructMatches(phospho_glycerol_pattern)
    if not phospho_matches:
        return False, "Phosphodiester group not properly attached to glycerol"

    # Ensure connectivity between phosphodiester group and cytidine moiety
    cytidine_matches = mol.GetSubstructMatches(cytidine_pattern)
    if not cytidine_matches:
        return False, "Cytidine moiety not properly attached to phosphodiester group"

    return True, "Molecule is a CDP-diacylglycerol"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17963',
        'name': 'CDP-diacylglycerol',
        'definition': 'A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.',
        'parents': []
    },
}