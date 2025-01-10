"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for detecting a steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected"

    # SMARTS pattern for detecting alpha 3-hydroxy group (stereochemistry aware)
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([OH])C")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "No 3alpha-hydroxy group detected"
    
    # Position 3 verification
    matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    three_alpha_hydroxyl = False
    for match in matches:
        # Check if the hydroxy group is at the 3-position of the steroid
        for atom_idx in match:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8:  # Oxygen atom
                carbon_idx = mol.GetAtomWithIdx(atom_idx).GetNeighbors()[0].GetIdx()
                if mol.GetAtomWithIdx(carbon_idx).GetSymbol() == 'C':
                    neighbors = mol.GetAtomWithIdx(carbon_idx).GetNeighbors()
                    if any(n.GetSymbol() == 'C' and n.GetIdx() != atom_idx for n in neighbors):
                        for n in neighbors:
                            if n.GetSymbol() == 'C':
                                two_bond_away = set(atom.GetIdx() for atom in n.GetNeighbors() if atom.GetSymbol() == 'C')
                                steroid_carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
                                if len(set(steroid_carbon_idxs).intersection(two_bond_away)) >= 8:
                                    # Assuming this means we're around position 3 in a steroid core
                                    three_alpha_hydroxyl = True
                                    break

    if not three_alpha_hydroxyl:
        return False, "Identified alpha hydroxy group is not at the 3-position"

    return True, "3alpha-hydroxy steroid structure identified"

# Example usage
smiles_example = "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)CO"
result, reason = is_3alpha_hydroxy_steroid(smiles_example)
print(f"Is 3alpha-hydroxy steroid: {result}, Reason: {reason}")