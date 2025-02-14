"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glyceride in which any one of the R groups (position not specified)
    is an acyl group attached via an ester linkage to glycerol, while the remaining two R groups
    can be either H or alkyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern (allowing for substitutions and stereochemistry)
    glycerol_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Define ester linkage pattern (generalized for any attachment point)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check that the ester is connected to the glycerol backbone
    ester_atom_idx = ester_matches[0][2]  # Index of the [C] in "C(=O)O[C]"
    connected_to_glycerol = False
    for match in glycerol_matches:
        if ester_atom_idx in match:
            connected_to_glycerol = True
            glycerol_atoms = match
            break
    if not connected_to_glycerol:
        return False, "Ester group not connected to glycerol backbone"

    # Check for additional ester or amide groups (should not have any)
    other_ester_pattern = Chem.MolFromSmarts("C(=O)[O,N][!#1]")
    total_ester_matches = mol.GetSubstructMatches(other_ester_pattern)
    if len(total_ester_matches) != 1:
        return False, f"Total ester/amide groups found: {len(total_ester_matches)}, expected exactly 1"

    # Exclude molecules with phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)[O]")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, not a monoacylglycerol"

    # Check the substituents on the glycerol backbone
    substituent_count = 0
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in glycerol_atoms:
                if neighbor.GetAtomicNum() != 1:  # Exclude hydrogen
                    substituent_count += 1
    # Since one position is occupied by acyl group, remaining can have up to 2 substituents
    if substituent_count > 2:
        return False, f"Too many substituents on glycerol backbone: found {substituent_count}, expected up to 2"

    return True, "Contains glycerol backbone with one acyl group attached via ester linkage"