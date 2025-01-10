"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are monoacylglycerol phosphates obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the glycerol phosphate backbone
    # This pattern matches glycerol backbone connected to a phosphate group
    glycerol_phosphate_smarts = "[C;X4][C;X4][C;X4]OP(=O)(O)O"
    glycerol_phosphate_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)

    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol phosphate backbone found"

    # Identify the glycerol backbone carbons
    glycerol_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    glycerol_carbons = set()
    for match in glycerol_matches:
        # Indices of the three glycerol carbons in the match
        glycerol_c_idx = match[0:3]
        glycerol_carbons.update(glycerol_c_idx)

    # Define a pattern for ester bonds (acyl chains attached via ester linkage)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;!$(C=O)]")  # Ester linkage not adjacent to another carbonyl
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Count ester bonds connected to glycerol carbons
    acyl_ester_count = 0
    for match in ester_matches:
        ester_oxygen_idx = match[2]
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        # Check if the ester oxygen is connected to one of the glycerol carbons
        for neighbor in ester_oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() in glycerol_carbons:
                acyl_ester_count += 1
                break

    if acyl_ester_count != 1:
        return False, f"Expected 1 acyl chain attached via ester bond to glycerol backbone, found {acyl_ester_count}"

    # Ensure there are no additional ester bonds (excluding phosphate ester bonds)
    total_acyl_esters = acyl_ester_count
    if total_acyl_esters > 1:
        return False, f"Found {total_acyl_esters} acyl ester bonds, expected only 1"

    # Check that the phosphate group is not further substituted (e.g., with choline or ethanolamine)
    # Define a pattern for substituted phosphate groups
    substituted_phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)[!O]")
    if mol.HasSubstructMatch(substituted_phosphate_pattern):
        return False, "Phosphate group is substituted (not a lysophosphatidic acid)"

    return True, "Molecule is a lysophosphatidic acid with glycerol phosphate backbone and one acyl chain"