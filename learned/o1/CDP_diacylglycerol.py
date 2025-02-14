"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import Draw

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is a glycerol backbone with two fatty acid chains attached via ester bonds
    at positions 1 and 2, and a CDP (cytidine diphosphate) group attached via a phosphodiester bond
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

    # Define SMARTS patterns

    # Glycerol backbone (with or without stereochemistry)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO)(CO)")  # Glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Ester linkage (fatty acid chains attached via ester bonds)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")  # Ester linkage to chiral carbon
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Phosphodiester linkage to CDP group
    # Phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphodiester linkage found"

    # Cytidine moiety pattern (includes ribose sugar and cytosine base)
    cytidine_pattern = Chem.MolFromSmarts("N1C=CC(=NC1=O)N[C@H]2O[C@@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine moiety found"

    # Check connectivity between glycerol backbone and phosphate group
    # Find the glycerol carbon connected to the phosphate
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    connected = False
    for g_match in glycerol_matches:
        glycerol_atom_indices = g_match
        for p_match in phosphate_matches:
            phosphate_atom_indices = p_match
            # Check bonds between glycerol oxygen and phosphate phosphorus
            for idx in glycerol_atom_indices:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == "O":
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in phosphate_atom_indices and nbr.GetSymbol() == "P":
                            connected = True
                            break
    if not connected:
        return False, "Phosphate group not connected to glycerol backbone"

    # Check connectivity between phosphate group and cytidine moiety
    cytidine_matches = mol.GetSubstructMatches(cytidine_pattern)
    if not cytidine_matches:
        return False, "Cytidine moiety not properly attached"

    return True, "Molecule is a CDP-diacylglycerol"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17963',
        'name': 'CDP-diacylglycerol',
        'definition': 'A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.',
        'parents': []
    },
}