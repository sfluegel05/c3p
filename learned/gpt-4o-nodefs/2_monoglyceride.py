"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride consists of a glycerol backbone, with a single fatty acid attached
    at the 2 position (middle) via an ester bond, and hydroxyl groups on the 1 and 3 positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the glycerol backbone with ester linkage at 2-position.
    glycerol_2_mono_pattern = Chem.MolFromSmarts("C(CO)(CO)OC(=O)[C]")
    if not mol.HasSubstructMatch(glycerol_2_mono_pattern):
        return False, "Structure does not match glycerol backbone with ester at 2-position pattern"

    # Check that the ester-part is a long chain (fatty acid like)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)OC(C)"))
    
    for match in ester_matches:
        # Focus on the carbon chain length after the ester linkage
        ester_carbon_idx = match[4]  # The atom in the ester
        hydrocarbon_chain_length = 0
        atom_queue = [mol.GetAtomWithIdx(ester_carbon_idx)]
        visited_atoms = set()
        
        while atom_queue:
            current_atom = atom_queue.pop(0)
            visited_atoms.add(current_atom.GetIdx())
            if current_atom.GetSymbol() == 'C' and current_atom.GetIdx() != ester_carbon_idx:
                hydrocarbon_chain_length += 1
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms:
                    atom_queue.append(neighbor)

        if hydrocarbon_chain_length >= 10:  # Demand at least 11 carbons
            return True, "Structure matches 2-monoglyceride with ester linkage at glycerol 2 position"
    
    return False, "No sufficient fatty acid-like esterification at glycerol 2 position found"