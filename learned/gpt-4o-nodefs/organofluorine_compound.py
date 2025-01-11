"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound contains one or more carbon-fluorine (C-F) bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over atoms to find C-F bonds
    has_c_f_bond = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 9:  # Check if atom is Fluorine (F)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Check if neighbor is Carbon (C)
                    has_c_f_bond = True
                    break
        if has_c_f_bond:
            break

    if has_c_f_bond:
        return True, "Contains one or more carbon-fluorine (C-F) bonds"

    return False, "No carbon-fluorine (C-F) bonds found"

# Example usage
example_smiles = "C[C@@H](OC1=NC(F)=C(Cl)C(N)=C1Cl)C(O)=O"  # Fluchloraminopyr
result, reason = is_organofluorine_compound(example_smiles)
print(result, reason)