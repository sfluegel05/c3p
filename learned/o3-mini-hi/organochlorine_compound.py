"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is defined as one that has a main (largest) organic fragment
            containing at least one carbon–chlorine bond.
Improvements:
  - Instead of simply looking for any C–Cl bond in the entire molecule, the code splits the molecule
    into fragments and uses the largest one.
  - The fragment must consist solely of atoms typically found in organic molecules.
  - An extra heuristic ensures that the fragment is “organic enough” by requiring a minimum fraction of carbon atoms
    (and a maximum ratio of non‐carbon to carbon heavy atoms) so that inorganic salts or very complex (e.g. peptide) structures
    are not mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is a typical organochlorine compound based on its SMILES.
    A typical organochlorine compound (by our improved definition) is one whose largest fragment is:
      - A valid molecule composed solely of atoms from a typical organic set.
      - Contains at least one carbon.
      - Contains at least one carbon–chlorine bond.
      - Has a sufficiently high fraction of carbon atoms (heuristic for typical organic compounds).
      - Does not have an excessive ratio of heteroatoms (non‐carbon heavy atoms) to carbon atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organochlorine compound, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split the molecule into fragments and work on the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found"
    largest_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Check that the largest fragment contains at least one carbon.
    if not any(atom.GetAtomicNum() == 6 for atom in largest_frag.GetAtoms()):
        return False, "Largest fragment does not contain any carbon atoms"
    
    # Define allowed elements (atomic numbers) typically found in organic molecules.
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}  # H, B, C, N, O, F, Si, P, S, Cl, Br, I
    for atom in largest_frag.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non‐organic element: {atom.GetSymbol()}"
    
    # Heuristic: ensure the fragment is “organic enough” by checking carbon enrichment.
    heavy_atoms = [atom for atom in largest_frag.GetAtoms() if atom.GetAtomicNum() > 1]
    num_heavy = len(heavy_atoms)
    num_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    if num_heavy == 0:
        return False, "No heavy atoms found"
    carbon_fraction = num_carbon / num_heavy
    # Require at least 40% of heavy atoms to be carbon.
    if carbon_fraction < 0.40:
        return False, f"Organic fraction too low (carbon fraction = {carbon_fraction:.2f})"
    
    # Extra heuristic: avoid fragments that have an unusually high ratio
    # of non-carbon atoms to carbon atoms. (Threshold can be tuned.)
    hetero_ratio = (num_heavy - num_carbon) / num_carbon
    if hetero_ratio > 2.5:
        return False, f"Excessive heteroatom content (non-C/C ratio = {hetero_ratio:.2f})"
    
    # Use a SMARTS pattern that matches a C–Cl bond.
    ccl_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
    if not largest_frag.HasSubstructMatch(ccl_pattern):
        return False, "Does not contain any carbon–chlorine bonds"
    
    return True, "Contains at least one carbon–chlorine bond and the main fragment is organic"

# Example usage:
if __name__ == "__main__":
    # Test a few SMILES:
    test_smiles = [
        "ClCCBr",
        "C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)CN3C=NN=N3)O[C@H]1CN(C)CC4=CC(=C(C=C4)Cl)Cl)[C@@H](C)CO",
        "ClC(Cl)(Cl)C(c1ccccc1)c1ccccc1",
        "C(C(CO)O)N1C=C(I)C(C(=C1)I)=O",  # iopydol; no C-Cl bond.
    ]
    for smi in test_smiles:
        result, reason = is_organochlorine_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*40}")