"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is defined as one that contains at least one carbon-chlorine bond.
Note: This improved version not only looks for a C–Cl bond (using a SMARTS pattern that works for both aliphatic and aromatic bonds)
      but also verifies that the molecule is “organic” by checking that it only contains atoms typically found in organic compounds.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound should:
      - Be a valid molecule.
      - Contain at least one carbon-chlorine bond.
      - Consist only of atoms typically found in organic compounds
        (for example: H, B, C, N, O, F, Si, P, S, Cl, Br, I).
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organochlorine compound, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a set of allowed atomic numbers for an organic molecule:
    # H (1), B (5), C (6), N (7), O (8), F (9), Si (14), P (15), S (16), Cl (17), Br (35), I (53)
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non‐organic element: {atom.GetSymbol()}"
    
    # Use an improved SMARTS that matches any bond (regardless of bond type) between a carbon (atomic number 6)
    # and chlorine (atomic number 17). The SMARTS "[#6][Cl]" will match both aliphatic and aromatic C-Cl bonds.
    ccl_pattern = Chem.MolFromSmarts("[#6][Cl]")
    if not mol.HasSubstructMatch(ccl_pattern):
        return False, "Does not contain any carbon-chlorine bonds"
    
    # (Optional additional filter:)
    # In some datasets many false positives were large natural products or ion complexes that—even though they contain a C–Cl bond—are not
    # considered typical organochlorine compounds. One might decide to check if the molecule’s size is within a typical range.
    # For example, if desired, uncomment the following lines to require a molecular weight less than 1000 Da.
    #
    # mol_wt = Descriptors.ExactMolWt(mol)
    # if mol_wt > 1000:
    #     return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a typical organochlorine compound"
    
    # If we made it here, the molecule is valid, organic, and has at least one C–Cl bond.
    return True, "Contains at least one carbon-chlorine bond and is an organic molecule"
    
# Example usage:
# test_smiles = "ClCCBr"  # 1-bromo-2-chloroethane example
# print(is_organochlorine_compound(test_smiles))