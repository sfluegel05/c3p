"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one chlorine atom within an organic molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    has_carbon = False
    has_chlorine = False

    # Iterate over all atoms in the molecule to check for carbon and chlorine
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:  # Carbon
            has_carbon = True
        elif atomic_num == 17:  # Chlorine
            has_chlorine = True

    # Check if both carbon and chlorine are present
    if has_carbon and has_chlorine:
        return True, "Contains both carbon and chlorine atoms"
    else:
        reason_parts = []
        if not has_carbon:
            reason_parts.append("No carbon atoms found")
        if not has_chlorine:
            reason_parts.append("No chlorine atoms found")
        return False, "; ".join(reason_parts)