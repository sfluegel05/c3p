"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is defined as one that contains at least one carbon-chlorine bond.
Improvement: Instead of rejecting the entire molecule if any atom lies outside a limited allowed set,
we split the molecule into fragments and apply the organic filter and substructure search to the largest fragment.
This helps avoid rejecting molecules that are salts or include additional ions.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound should:
      - Be a valid molecule.
      - Have a main (largest) fragment that is organic (i.e. composed only of typically organic atoms).
      - Contain at least one carbon-chlorine bond in that fragment.
    
    The allowed organic elements include: H, B, C, N, O, F, Si, P, S, Cl, Br, I
    (atomic numbers: 1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is considered an organochlorine compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For molecules with multiple fragments, analyze only the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found"
    # Choose the fragment with the maximum number of heavy atoms.
    largest_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Check that the largest fragment contains at least one carbon atom.
    if not any(atom.GetAtomicNum() == 6 for atom in largest_frag.GetAtoms()):
        return False, "Largest fragment does not contain any carbon atoms"
    
    # Define allowed atomic numbers for an organic fragment.
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    for atom in largest_frag.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non‐organic element: {atom.GetSymbol()}"
    
    # Use a SMARTS pattern that matches a bond (regardless of bonding type) between a carbon and chlorine.
    # The pattern "[#6]-[Cl]" ensures an explicit bond between a carbon (atomic number 6) and chlorine (atomic number 17).
    ccl_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
    if not largest_frag.HasSubstructMatch(ccl_pattern):
        return False, "Does not contain any carbon-chlorine bonds"
    
    # (Optional filter: Only consider molecules with a typical size for organochlorine compounds.)
    # mol_wt = Descriptors.ExactMolWt(largest_frag)
    # if mol_wt > 1000:
    #     return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a typical organochlorine compound"
    
    return True, "Contains at least one carbon-chlorine bond and the main fragment is organic"

# Example usage:
# test_smiles = "ClCCBr"  # 1-bromo-2-chloroethane example
# print(is_organochlorine_compound(test_smiles))