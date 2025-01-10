"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are monoacylglycerol phosphates obtained by hydrolytic removal 
    of one of the two acyl groups of any phosphatidic acid or derivatives therein.

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

    # Ensure molecule has a phosphate group connected to glycerol
    glycerol_phosphate_smarts = "[C@@H](CO[P](=O)(O)O)O"
    gp_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(gp_pattern):
        return False, "No glycerol phosphate backbone found"

    # Check for absence of substituents on the phosphate group (no additional atoms attached to phosphate)
    phosphate_atom = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            phosphate_atom = atom
            break
    if phosphate_atom is None:
        return False, "No phosphate group found"

    phosphate_neighbors = phosphate_atom.GetNeighbors()
    num_non_oxygen_neighbors = sum(1 for neighbor in phosphate_neighbors if neighbor.GetAtomicNum() != 8)
    if num_non_oxygen_neighbors > 1:
        return False, "Phosphate group is substituted (not a lysophosphatidic acid)"

    # Identify glycerol backbone carbons
    glycerol_carbons = []
    for match in mol.GetSubstructMatches(gp_pattern):
        glycerol_carbons.extend([match[0], match[1], match[2]])  # Indices of the glycerol carbons
        break  # Assuming one glycerol phosphate backbone

    # Define ester bond pattern (acyl chain attached via ester linkage to glycerol)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;H2][C;H][C;H2]")  # Acyl ester attached to glycerol carbons
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Check the number of acyl chains attached via ester bonds to glycerol carbons
    acyl_ester_count = 0
    for match in ester_matches:
        ester_oxygen_idx = match[2]
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        # Check if the ester oxygen is connected to one of the glycerol carbons
        for neighbor in ester_oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() in glycerol_carbons:
                acyl_ester_count += 1
                break

    if acyl_ester_count == 0:
        return False, "No acyl chains attached via ester bonds to glycerol backbone"
    elif acyl_ester_count > 2:
        return False, f"Expected 1 or 2 acyl chains attached via ester bonds to glycerol backbone, found {acyl_ester_count}"

    # Check for additional acyl chains attached via ether bonds (LPAs should not have ether-linked chains)
    ether_pattern = Chem.MolFromSmarts("[C;H2][O][C;H2][C;H][C;H2]")  # Ether linkage to glycerol carbons
    if mol.HasSubstructMatch(ether_pattern):
        return False, "Found ether-linked chains, which are not typical for LPAs"

    return True, "Molecule is a lysophosphatidic acid with glycerol phosphate backbone and acyl chain(s)"