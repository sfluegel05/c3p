"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Corrected SMARTS pattern for sulfonic acid group connected via carbon-sulfur bond
    # Sulfonic acid group: S(=O)(=O)-O[-] attached to carbon via S-C bond
    sulfonic_acid_pattern = Chem.MolFromSmarts("[C]-S(=O)(=O)-[O;H1,-1]")
    if sulfonic_acid_pattern is None:
        return None, "Invalid sulfonic acid SMARTS pattern"

    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group connected via carbon-sulfur bond found"

    # Check for lipid characteristics
    # Instead of using unsupported SMARTS quantifiers, programmatically find long carbon chains
    def get_longest_carbon_chain(mol):
        max_length = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                visited = set()
                stack = [(atom, 1)]
                while stack:
                    current_atom, length = stack.pop()
                    idx = current_atom.GetIdx()
                    if idx in visited:
                        continue
                    visited.add(idx)
                    if length > max_length:
                        max_length = length
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                            stack.append((neighbor, length + 1))
        return max_length

    longest_chain = get_longest_carbon_chain(mol)
    if longest_chain < 12:
        return False, f"No long hydrocarbon chain of at least 12 carbons found (longest chain: {longest_chain})"

    return True, "Contains sulfonic acid group connected via carbon-sulfur bond and lipid characteristics"