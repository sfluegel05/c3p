"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (a C25 precursor), including compounds
where the C25 skeleton has been rearranged or partially truncated (often via the loss of methyl groups).
Note: This program uses multiple heuristics (scaffold carbon count, aromatic ring count,
and full-molecule carbon count plus aromatic content) to reduce both false positives and false negatives.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Primary decision is made via a heuristic:
      1. Compute the Bemis-Murcko scaffold and count its carbon atoms.
         Terpenoids from a C25 precursor (5 isoprene units) are expected to yield a scaffold
         with roughly 22 to 28 carbons. In addition, most terpenoid frameworks are largely aliphatic.
         Therefore, we also count the number of fully aromatic rings in the scaffold.
      2. If the scaffold does not land in that range – for example, due to extensive modifications –
         we fall back on computing the total number of carbons in the full molecule and assessing its
         aromatic content.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a sesterterpenoid, False otherwise.
        str: Explanation of the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Compute the Bemis–Murcko scaffold
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating scaffold: {e}"
    
    if scaffold is None:
        return False, "Could not compute Murcko scaffold for molecule."
    
    # Count carbon atoms in scaffold
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count fully aromatic rings in scaffold.
    ring_info = scaffold.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        # Check if every atom in the ring is aromatic.
        if all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Primary heuristic: if the scaffold has a roughly C22-28 backbone and is mainly aliphatic.
    if 22 <= scaffold_carbons <= 28 and aromatic_ring_count <= 1:
        return True, f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s), consistent with a sesterterpenoid backbone."
    
    # Otherwise, as fallback, examine the full molecule.
    full_mol_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if full_mol_carbons == 0:
        return False, "Full molecule has no carbon atoms."
    
    # Count aromatic carbons in the full molecule.
    aromatic_carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    aromatic_fraction = aromatic_carbon_count / full_mol_carbons

    # Accept as a sesterterpenoid if the full molecule has a carbon count near 25 (allowing a slight range)
    # and a low aromatic fraction (here, less than 30% aromatics).
    if 22 <= full_mol_carbons <= 30 and aromatic_fraction < 0.3:
        return True, f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic carbons, consistent with a modified sesterterpenoid skeleton."
    
    return False, f"Scaffold has {scaffold_carbons} carbons (and {aromatic_ring_count} aromatic ring(s)) or full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic content, which is not consistent with a typical sesterterpenoid."

# Example usage:
if __name__ == "__main__":
    # (2Z,6E,10E,14E)-geranylfarnesol is a classic example of a sesterterpenoid precursor (C25)
    test_smiles = "C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)/C(=C\\CO)/C"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)