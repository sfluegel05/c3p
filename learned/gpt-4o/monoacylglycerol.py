"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone and one acyl group (ester bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced glycerol backbone pattern allowing for ester derivation
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Ester linkage: -O-C(=O)-
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 1"
    
    # Ensure the acyl group is sufficiently long (e.g., 4+ carbons for this classification)
    acyl_start = ester_matches[0][1]  # Carbonyl carbon position
    carbon_count = 0
    visited_atoms = set()

    def count_acyl_carbons(start_atom):
        nonlocal carbon_count
        current_atom = start_atom
        while current_atom and current_atom.GetIdx() not in visited_atoms:
            if current_atom.GetAtomicNum() == 6:  # Check for carbon
                carbon_count += 1
                visited_atoms.add(current_atom.GetIdx())
                # Move to the next carbon in chain, skipping ester oxygen
                next_atoms = [atom for atom in current_atom.GetNeighbors()
                              if atom.GetAtomicNum() == 6 and atom.GetIdx() not in visited_atoms]
                current_atom = next_atoms[0] if next_atoms else None
            else:
                break

    # Find first carbon in acyl chain
    acyl_neighs = [atom for atom in mol.GetAtomWithIdx(acyl_start).GetNeighbors() if atom.GetAtomicNum() == 6]
    if acyl_neighs:
        count_acyl_carbons(acyl_neighs[0])

    # Check carbon count in acyl chain
    if carbon_count < 4:
        return False, f"Acyl chain length is {carbon_count}, which is too short for typical monoacylglycerol"

    return True, "Contains glycerol backbone with one acyl group attached via ester bond"