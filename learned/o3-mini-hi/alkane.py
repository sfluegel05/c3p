"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: Alkane
Definition: An acyclic branched or unbranched hydrocarbon having the general formula CnH2n+2,
and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic saturated hydrocarbon with the general formula CnH2n+2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to reliably count them
    mol = Chem.AddHs(mol)

    # Check that the molecule contains only carbon and hydrogen atoms.
    allowed_atomic_nums = {1, 6}  # Hydrogen and Carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Atom '{atom.GetSymbol()}' found which is not carbon or hydrogen"

    # Count the number of Carbon and Hydrogen atoms.
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_hydrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    # Check if formula matches CnH2n+2.  (For n==0, not an alkane)
    if n_carbons == 0:
        return False, "No carbon atoms found"
    if n_hydrogens != 2 * n_carbons + 2:
        return False, f"Hydrogen count does not match formula CnH2n+2 (found C{n_carbons}H{n_hydrogens})"

    # Check that all bonds in the molecule are single bonds ensuring saturation.
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False, "Non-single bond found, molecule is not fully saturated"

    # Check that the molecule is acyclic (no rings).
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings which is not allowed in alkanes"

    return True, "Molecule is an acyclic, saturated hydrocarbon with formula CnH2n+2"

# Example usage (uncomment to test):
# smiles_examples = [
#     "C(C(CC)C)CCCCCCCCC",  # 3-methyl-Tridecane
#     "[H]C([H])([H])[H]",  # methane
#     "CC(C)C",             # isobutane
#     "CCCC"                # butane
# ]
# for s in smiles_examples:
#     result, reason = is_alkane(s)
#     print(f"SMILES: {s} => {result} ({reason})")