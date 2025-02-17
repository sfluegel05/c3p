"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHALCOGEN - Any p-block element belonging to the group 16 family of the periodic table.
Valid chalcogens: Oxygen (O), Sulfur (S), Selenium (Se), Tellurium (Te), Polonium (Po)
A valid chalcogen entry should be a single atom as represented by its SMILES string.
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen element based on its SMILES string.
    A chalcogen is defined as any p-block element belonging to the group 16 family of the periodic table,
    i.e., oxygen, sulfur, selenium, tellurium, and polonium.
    Additionally, to be classified as a pure chalcogen, the SMILES should represent a single isolated atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a single-atom chalcogen, False otherwise
        str: Reason for classification
    """
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule consists of exactly one atom.
    # Note: Explicit hydrogens are included in GetNumAtoms(), so we also ensure there are no extra atoms.
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not represent a single isolated atom"
    
    # Get the one atom in the molecule.
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    
    # Define allowed atomic numbers for group 16 chalcogens:
    # Oxygen (8), Sulfur (16), Selenium (34), Tellurium (52), Polonium (84)
    allowed_atomic_nums = {8, 16, 34, 52, 84}
    if atomic_num in allowed_atomic_nums:
        return True, f"Element with atomic number {atomic_num} is a chalcogen"
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen"

# Example test cases (uncomment for testing):
# test_smiles = [
#     "[218Po]", "[202Po]", "[208Po]", "[209Po]", "[190Po]", "[201Po]",
#     "[19O]", "[216Po]", "[38S]", "[36S]", "[203Po]", "[Te]",
#     "[210Po]", "[18O]", "[Po]", "[217Po]", "[82Se]", "[Se]",
#     "[206Po]", "[211Po]", "[193Po]", "[194Po]", "[195Po]", "[197Po]",
#     "[33S]", "[125Te]", "[204Po]", "[S]", "[200Po]", "[15O]",
#     "[196Po]", "[35S]", "[37S]", "[199Po]", "[191Po]", "[198Po]",
#     "[34S]", "[17O]", "[214Po]",
#     # False positive examples:
#     "[O+][H]", "[S--]", "[H][S+][H]", "[H]O[H]", "[H][Se][H]",
#     "[Te+][H]", "[Se--]", "[Se-]", "[H][Te-][H]", "[Se+][H]",
#     "[S-2]", "[S-]", "[H][Te][H]", "[Se++]", "[H]S([H])([H])[H]",
#     "[H][Te+]([H])[H]", "[O-][H]", "[Te--]", "[S+2]", "[H][Po][H]",
#     "[H][S+]([H])[H]", "[S+][H]", "[Se][H]", "[H][S-]([H])[H]", "[O-]"
# ]
# for ts in test_smiles:
#     status, reason = is_chalcogen(ts)
#     print(f"{ts}: {status} | {reason}")