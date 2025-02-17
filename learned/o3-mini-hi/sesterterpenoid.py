"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (a C25 precursor) including compounds 
where the C25 skeleton has been rearranged or partly truncated.
This script uses a combination of Bemis–Murcko scaffold analysis and full‐molecule statistics 
(including carbon count, aromaticity, sp3 fraction, and the ratio of scaffold to full molecule carbons)
to decide whether a compound is a sesterterpenoid.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    
    Improved Heuristics:
      • Compute the Bemis–Murcko scaffold.
      • Count carbon atoms in the scaffold and in the full molecule.
      • Count fully aromatic rings in the scaffold and overall aromatic fraction.
      • Compute the fraction of sp3-hybridized carbons.
      • Compute the ratio of scaffold carbons to full-molecule carbons.
      
    We expect a true sesterterpenoid (or a modified derivative) to have:
      - A terpenoid-like scaffold between roughly 18 and 32 carbons.
      - A scaffold that makes up at least 65% of the full molecule’s carbon count.
      - Very low overall aromatic content (<30% aromatic carbons).
      - High sp3 character (sp3 fraction ≥0.6 when looking at the scaffold or ≥0.5 in the full molecule).
      
    We allow a fallback (full molecule-based) if the scaffold conditions are not fully met 
    but overall the full molecule falls in reasonable windows (22–60 carbons, low aromaticity, decent sp3)
    and the scaffold still represents a significant fraction.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a sesterterpenoid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Try to compute the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating scaffold: {e}"
    if scaffold is None:
        return False, "Could not compute Murcko scaffold for molecule."
    
    # Count carbons in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count fully aromatic rings in the scaffold.
    aromatic_ring_count = 0
    ring_info = scaffold.GetRingInfo()
    for ring in ring_info.AtomRings():
        if all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Count full molecule carbons and sp3 carbons.
    full_atoms = mol.GetAtoms()
    full_mol_carbons = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6)
    sp3_count = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6 and 
                     atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    sp3_fraction = sp3_count / full_mol_carbons if full_mol_carbons > 0 else 0

    # Count aromatic carbons in the full molecule.
    aromatic_carbon_count = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    aromatic_fraction = aromatic_carbon_count / full_mol_carbons if full_mol_carbons > 0 else 0
    
    # Compute ratio: what fraction of the full carbon count is in the scaffold.
    scaffold_ratio = scaffold_carbons / full_mol_carbons if full_mol_carbons > 0 else 0

    # Primary heuristic based on the scaffold.
    # We expect a modified sesterterpenoid to have an underlying terpene core:
    # - Scaffold carbons between 18 and 32,
    # - Minimal aromatic rings (≤1),
    # - High sp3 fraction (≥0.6),
    # - And the scaffold should make up at least 65% of full molecule carbons.
    if (18 <= scaffold_carbons <= 32 and aromatic_ring_count <= 1 and sp3_fraction >= 0.6
        and scaffold_ratio >= 0.65):
        return True, (f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s), "
                      f"sp3 fraction of {sp3_fraction:.0%}, and represents {scaffold_ratio:.0%} of the full "
                      f"{full_mol_carbons} carbons, which is consistent with a sesterterpenoid backbone.")
    
    # Fallback heuristic: use full molecule statistics with a wider window.
    # Here we allow full molecule carbon counts roughly between 22 and 60,
    # low aromatic carbon fraction (<30%),
    # moderate sp3 fraction (≥0.5),
    # and require that the scaffold represents at least 65% of the total carbons.
    if (22 <= full_mol_carbons <= 60 and aromatic_fraction < 0.3 
        and sp3_fraction >= 0.5 and scaffold_ratio >= 0.65):
        return True, (f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} "
                      f"aromatic carbon content and sp3 fraction of {sp3_fraction:.0%}; "
                      f"scaffold ({scaffold_carbons} carbons) represents {scaffold_ratio:.0%} of all carbons. "
                      f"These characteristics are consistent with a modified sesterterpenoid skeleton.")
    
    # If neither primary nor fallback conditions are met, the compound is not classified as sesterterpenoid.
    return False, (f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s) "
                   f"and makes up {scaffold_ratio:.0%} of {full_mol_carbons} total carbons; "
                   f"full molecule has {aromatic_fraction:.0%} aromatic content and sp3 fraction of {sp3_fraction:.0%}. "
                   f"These values are not consistent with a typical sesterterpenoid.")

# Example usage:
if __name__ == "__main__":
    # Test the function with one of the provided sesterterpenoid examples.
    test_smiles = "O=C1[C@@H]2[C@@H](C=O)CC[C@H]3[C@@](O)([C@H](CCCC(C)C)C)CC[C@@]3(C[C@@H]2[C@@](C1)(O)C)C"  # Zizanin-B
    result, reason = is_sesterterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)