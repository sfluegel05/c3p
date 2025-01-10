"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid
    with a chain length of less than C6 and no non-hydrocarbon substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES, return False if invalid
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid group SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for aromatic rings using a workaround method by checking bond orders
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            return False, "Contains aromatic rings, not aliphatic"

    # Ensure molecule only contains C, H, O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:  # Hydrogen, Carbon, Oxygen
            return False, "Contains non-hydrocarbon substituents"

    # Verify the longest carbon chain connected to carboxylic acid
    for match in mol.GetSubstructMatches(carboxylic_acid_pattern):
        carboxylic_c_index = match[0]  # Carbon atom of the carboxyl group
        visited = set()
        to_visit = [(carboxylic_c_index, 0)]  # Start with the carboxyl carbon, depth 0
        
        max_chain_length = 0
        while to_visit:
            atom_idx, chain_length = to_visit.pop(0)
            max_chain_length = max(max_chain_length, chain_length)
            visited.add(atom_idx)
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if neighbor.GetIdx() not in visited:  # Avoid cycles, prevents revisiting
                    # Track only carbon atoms for chain length
                    if neighbor.GetAtomicNum() == 6:
                        to_visit.append((neighbor.GetIdx(), chain_length + 1))

        # Maximum permissible carbon chain length is less than 6 (excluding the carboxyl carbon itself)
        if max_chain_length >= 5:
            return False, f"Carbon chain too long, found {max_chain_length} carbons"

    # Check the molecular weight to ensure it stays within expected range (optional sanity check)
    if MolWt(mol) > 150:
        return False, "Molecular weight too high for a short-chain fatty acid"

    return True, "Valid short-chain fatty acid with appropriate carbon chain length"

__metadata__ = {
    'chemical_class': {
        'name': 'short-chain fatty acid',
        'definition': 'An aliphatic monocarboxylic acid with a chain length of less than C6. If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.'
    }
}