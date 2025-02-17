"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHALCOGEN - Any p-block element belonging to the group 16 family of the periodic table.
Valid chalcogens: Oxygen (O), Sulfur (S), Selenium (Se), Tellurium (Te), Polonium (Po)
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen element based on its SMILES string.
    A chalcogen is defined as any p-block element belonging to the group 16 family of the periodic table,
    i.e., oxygen, sulfur, selenium, tellurium, and polonium.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a chalcogen atom, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that the molecule contains exactly one atom.
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not represent a single atom"
    
    # Obtain the atom
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    
    # Define allowed atomic numbers for group 16: O (8), S (16), Se (34), Te (52), Po (84)
    allowed_atomic_nums = {8, 16, 34, 52, 84}
    
    if atomic_num in allowed_atomic_nums:
        return True, f"Element with atomic number {atomic_num} is a chalcogen"
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen"

# Example test cases (uncomment for testing):
# test_smiles = ["[218Po]", "[19O]", "[38S]", "[Se]", "[125Te]", "[Fe]"]
# for ts in test_smiles:
#     status, reason = is_chalcogen(ts)
#     print(f"{ts}: {status} | {reason}")