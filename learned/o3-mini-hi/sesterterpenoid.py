"""
Classifies: CHEBI:26660 sesterterpenoid
"""
# sesterterpenoid.py
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (a C25 precursor), including compounds
where the C25 skeleton has been rearranged or partially truncated (often via the loss of methyl groups).
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    The decision is made by computing the molecule's Bemis–Murcko scaffold and counting its
    number of carbon atoms. Sesterterpenoids are derived from a C25 precursor (5x isoprene units)
    although slight modifications (eg removal of methyl groups) are possible. Here we use a range of
    roughly 22 to 28 carbons on the scaffold as a heuristic.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a sesterterpenoid, False otherwise.
        str: Explanation of the result.
    """
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to get the Bemis–Murcko scaffold of the molecule,
    # which represents the core framework (rings and linkers)
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating scaffold: {e}"
    
    if scaffold is None:
        return False, "Could not compute Murcko scaffold for molecule"

    # Count the number of carbon atoms in the scaffold.
    carbons_in_scaffold = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Use heuristic: if the scaffold has roughly between 22 and 28 carbons, 
    # then it is consistent with a C25 backbone (allowing for minor modifications).
    if 22 <= carbons_in_scaffold <= 28:
        return True, f"Scaffold has {carbons_in_scaffold} carbon atoms, consistent with a sesterterpenoid backbone."
    else:
        return False, f"Scaffold has {carbons_in_scaffold} carbon atoms, which is not consistent with a typical sesterterpenoid (expected near 25)."

# Example usage:
if __name__ == "__main__":
    # (2Z,6E,10E,14E)-geranylfarnesol is a classic example of a sesterterpenoid precursor (C25)
    test_smiles = "C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)/C(=C\\CO)/C"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)