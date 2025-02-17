"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is defined as a compound 
            whose main (largest) organic fragment contains at least one carbon–chlorine bond.
Improvements:
  - We split the molecule into fragments and pick the largest fragment.
  - We require that fragment is built from typical organic elements 
    and (if sufficiently large) that it is “organic‐enough” (by carbon fraction criteria).
  - We then search for a C–Cl bond in the largest fragment.
Note: For very small molecules (with few heavy atoms) the carbon fraction check is skipped.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES.
    For our purposes, the molecule is classified as an organochlorine compound if its largest
    fragment is comprised solely of typical organic atoms, contains at least one carbon,
    and in that fragment there is at least one C–Cl bond.
    
    In order to avoid excluding very small molecules (e.g. chloroform) we only enforce a minimum
    carbon fraction (and maximum heteroatom ratio) when the fragment has at least 10 heavy atoms.
      
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Split molecule into fragments and work with the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found"
    largest_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    # Check that the largest fragment contains at least one carbon.
    if not any(atom.GetAtomicNum() == 6 for atom in largest_frag.GetAtoms()):
        return False, "Largest fragment contains no carbon atoms"
    
    # Define allowed elements (atomic numbers) that are typically found in organic compounds.
    # (H, B, C, N, O, F, Si, P, S, Cl, Br, I)
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    for atom in largest_frag.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non‐organic element: {atom.GetSymbol()}"

    # Heuristic: Check that the fragment is "organic enough" when it is sufficiently large.
    heavy_atoms = [atom for atom in largest_frag.GetAtoms() if atom.GetAtomicNum() > 1]
    num_heavy = len(heavy_atoms)
    num_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)

    # Only enforce the following heuristics if the fragment is large.
    if num_heavy >= 10:
        carbon_fraction = num_carbon / num_heavy
        if carbon_fraction < 0.40:
            return False, f"Organic fraction too low (carbon fraction = {carbon_fraction:.2f})"
        hetero_ratio = (num_heavy - num_carbon) / num_carbon if num_carbon > 0 else float('inf')
        if hetero_ratio > 2.5:
            return False, f"Excessive heteroatom content (non-C/C ratio = {hetero_ratio:.2f})"

    # Use a SMARTS pattern to check for a C–Cl bond.
    ccl_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
    if not largest_frag.HasSubstructMatch(ccl_pattern):
        return False, "Does not contain any carbon–chlorine bonds"
    
    return True, "Contains at least one carbon–chlorine bond in the main organic fragment"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "ClCCBr",  # 1-bromo-2-chloroethane
        "C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)CN3C=NN=N3)O[C@H]1CN(C)CC4=CC(=C(C=C4)Cl)Cl)[C@@H](C)CO",  # a complex organochlorine
        "ClC(Cl)(Cl)C(c1ccccc1)c1ccccc1",  # 1,1,1-trichloro-2,2-diphenylethane
        "C(C(CO)O)N1C=C(I)C(C(=C1)I)=O",  # iopydol; does not contain any C-Cl bond
        "[H]C(Cl)(Cl)Cl",  # chloroform (small molecule; organic fraction check skipped)
        "FC(F)(F)[C@@H](Cl)Br",  # (S)-halothane
        "ClC(Br)Br",  # Chlorodibromomethane (valid if one accepts any C–Cl bond)
        "O=C(O)CC(Cl)C(O)=O",  # 2-chlorosuccinic acid
    ]
    for smi in test_smiles:
        result, reason = is_organochlorine_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*50}")