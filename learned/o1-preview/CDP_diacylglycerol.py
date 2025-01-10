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

    # Suppress sanitization warnings
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)

    # Define SMARTS patterns

    # Glycerol backbone: three carbons connected sequentially
    glycerol_pattern = Chem.MolFromSmarts("C(C)C")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Ester linkage pattern: O=C-O-C (ester group)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;!$(C=O)]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Phosphodiester linkage pattern: O-P(=O)(O)-O
    phosphodiester_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphodiester_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    if not phosphodiester_matches:
        return False, "No phosphodiester group found"

    # Cytidine moiety pattern: pyrimidine ring with amine and keto groups
    cytidine_pattern = Chem.MolFromSmarts("n1ccnc(N)c1")
    cytidine_matches = mol.GetSubstructMatches(cytidine_pattern)
    if not cytidine_matches:
        return False, "No cytidine base found"

    # Ribose sugar attached to cytidine base
    ribose_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)O1")
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "No ribose sugar found"

    # Check that the cytidine base is connected to ribose
    cytidine_atom_indices = [atom_idx for match in cytidine_matches for atom_idx in match]
    ribose_atom_indices = [atom_idx for match in ribose_matches for atom_idx in match]
    cytidine_neighbors = set()
    for idx in cytidine_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            cytidine_neighbors.add(neighbor.GetIdx())
    if not any(idx in ribose_atom_indices for idx in cytidine_neighbors):
        return False, "Cytidine base not connected to ribose sugar"

    # Check that the ribose is connected to the phosphodiester group
    phospho_atom_indices = [atom_idx for match in phosphodiester_matches for atom_idx in match]
    ribose_neighbors = set()
    for idx in ribose_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            ribose_neighbors.add(neighbor.GetIdx())
    if not any(idx in phospho_atom_indices for idx in ribose_neighbors):
        return False, "Ribose sugar not connected to phosphodiester group"

    # Check that the phosphodiester group is connected to the glycerol backbone
    glycerol_atom_indices = [atom_idx for match in glycerol_matches for atom_idx in match]
    phospho_neighbors = set()
    for idx in phospho_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            phospho_neighbors.add(neighbor.GetIdx())
    if not any(idx in glycerol_atom_indices for idx in phospho_neighbors):
        return False, "Phosphodiester group not connected to glycerol backbone"

    # Check that the glycerol backbone has two ester linkages
    ester_oxygen_indices = [match[2] for match in ester_matches]  # Oxygen connected to glycerol
    ester_glycerol_connections = 0
    for idx in ester_oxygen_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in glycerol_atom_indices:
                ester_glycerol_connections += 1
                break
    if ester_glycerol_connections < 2:
        return False, "Ester groups not properly attached to glycerol backbone"

    return True, "Molecule is a CDP-diacylglycerol"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17963',
        'name': 'CDP-diacylglycerol',
        'definition': 'A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.',
        'parents': []
    },
}