"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHALCOGEN - Any p-block element belonging to the group 16 family of the periodic table.
Only molecules that are represented as a single isolated atom (with no explicit additional atoms)
are accepted. Allowed chalcogens: Oxygen (atomic number 8), Sulfur (16), Selenium (34),
Tellurium (52), and Polonium (84).
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is represented as a neutral, single isolated atom that is
    a group 16 (chalcogen) element based on its SMILES string.
    
    The molecule must be written as a single atom (for example, "[Se]" rather than "[H][Se][H]").
    This ensures that species like water, selane, or tellane, where additional atoms are explicitly
    provided, are not misclassified.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a valid single-atom chalcogen, False otherwise.
        str: A reason describing the result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that the molecule is a single isolated atom.
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not represent a single isolated atom"
    
    # Get the sole atom.
    atom = mol.GetAtomWithIdx(0)
    
    # Check that the atom has no formal charge.
    if atom.GetFormalCharge() != 0:
        return False, f"Atom has a formal charge of {atom.GetFormalCharge()}, not neutral"
    
    # Allowed atomic numbers for group 16 (chalcogens).
    allowed_atomic_nums = {8, 16, 34, 52, 84}
    atomic_num = atom.GetAtomicNum()
    if atomic_num in allowed_atomic_nums:
        return True, f"Element with atomic number {atomic_num} is a chalcogen"
    else:
        return False, f"Element with atomic number {atomic_num} is not a chalcogen"

# Example usage (uncomment for testing):
# test_smiles = [
#     "[218Po]", "[202Po]", "[208Po]", "[209Po]", "[190Po]", "[201Po]",
#     "[19O]", "[216Po]", "[38S]", "[36S]", "[203Po]", "[Te]",
#     "[210Po]", "[18O]", "[Po]", "[217Po]", "[82Se]", "[Se]",
#     "[206Po]", "[211Po]", "[193Po]", "[194Po]", "[195Po]", "[197Po]",
#     "[33S]", "[125Te]", "[204Po]", "[S]", "[200Po]", "[15O]",
#     "[196Po]", "[35S]", "[37S]", "[199Po]", "[191Po]", "[198Po]",
#     "[34S]", "[17O]", "[214Po]",
#     # False positive examples: these should be rejected since they are not strictly single-atom.
#     "[H]O[H]", "[H][Se][H]", "[H][Te][H]", "[H]S([H])([H])[H]", "[H][Po][H]"
# ]
#
# for ts in test_smiles:
#     status, reason = is_chalcogen(ts)
#     print(f"{ts}: {status} | {reason}")