"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHALCOGEN - Any p-block element belonging to the group 16 family of the periodic table.
Valid chalcogens for this classification: Oxygen (O), Sulfur (S), Selenium (Se), Tellurium (Te), and Polonium (Po).
A valid entry must be a neutral, single-atom species.
"""

from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a neutral, single-atom chalcogen based on its SMILES string.
    
    A chalcogen is defined as any p-block element in group 16:
      - Oxygen (atomic number 8)
      - Sulfur (atomic number 16)
      - Selenium (atomic number 34)
      - Tellurium (atomic number 52)
      - Polonium (atomic number 84)
    
    The molecule must be exactly one atom and must have a formal charge of 0.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a neutral, single-atom chalcogen, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that no bonds are present (i.e., the molecule is a single isolated atom).
    # Note: Sometimes hydrogens are implicit; we rely on the fact that explicit atoms (even hydrogens)
    # would increase the atom count.
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not represent a single isolated atom"
    
    # Get the only atom in the molecule.
    atom = mol.GetAtomWithIdx(0)
    
    # Check formal charge; we require a neutral (zero charge) atom.
    if atom.GetFormalCharge() != 0:
        return False, f"Atom has a formal charge of {atom.GetFormalCharge()}, not neutral"
    
    # Define allowed atomic numbers for group 16 chalcogens.
    allowed_atomic_nums = {8, 16, 34, 52, 84}
    atomic_num = atom.GetAtomicNum()
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
#
# for ts in test_smiles:
#     status, reason = is_chalcogen(ts)
#     print(f"{ts}: {status} | {reason}")