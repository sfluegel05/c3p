"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is characterized by a glycerol backbone with two hydroxy groups acylated (esterified).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the glycerol backbone (modified with esters): OCC system
    glycerol_backbone_pattern = Chem.MolFromSmarts("OCC(O)C")
    matches = mol.GetSubstructMatches(glycerol_backbone_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Identify ester group: O=C-O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_connected_counts = []

    # Verify ester linkages are connected to the correct glycerol backbone carbons
    for match in matches:
        glycerol_carbon_indices = match[1:3]  # Central carbons in glycerol after oxygens
        esters_connected_to_backbone = 0

        for ester in ester_matches:
            ester_oxygen = ester[2]
            connected_atoms = mol.GetAtomWithIdx(ester_oxygen).GetNeighbors()
            # Try to find connection to one of the glycerol carbon atoms
            if any(atom.GetIdx() in glycerol_carbon_indices for atom in connected_atoms):
                esters_connected_to_backbone += 1

        ester_connected_counts.append(esters_connected_to_backbone)

    if 2 in ester_connected_counts:  # At least one match must have 2 connected esters
        return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds"

    return False, "Did not find two ester linkages properly connected to the glycerol backbone"