"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is classified as a wax based on its SMILES string.
    Waxes are characterized by long-chain fatty acids esterified to long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define ester bond pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    # Check for presence of ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond present"

    # Helper function to check for long aliphatic chains
    def is_long_aliphatic_chain(atom_idx):
        visited_atoms = set()
        stack = [atom_idx]

        carbon_count = 0
        while stack:
            atom_id = stack.pop()
            atom = mol.GetAtomWithIdx(atom_id)
            if atom.GetAtomicNum() == 6:  # Carbon
                carbon_count += 1
                
                if carbon_count >= 14:
                    return True

                visited_atoms.add(atom_id)
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited_atoms:
                        stack.append(neighbor_idx)

        return False

    # Iterate over ester bonds found
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for match in ester_matches:
        # Check the carbon chain lengths on both side of ester linkage
        carbonyl_carbon_idx = match[0]
        aliphatic_oxygen_idx = match[2]

        if is_long_aliphatic_chain(carbonyl_carbon_idx) and is_long_aliphatic_chain(aliphatic_oxygen_idx):
            return True, "Molecule contains long aliphatic ester typical of wax"

    return False, "Ester bond found but does not form a typical wax structure (long aliphatic chains missing)"