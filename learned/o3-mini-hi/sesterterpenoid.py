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
Improved heuristics below adjust thresholds (e.g. lowering the scaffold carbon count minimum from 18 to 14,
reducing the required scaffold-to-full carbon ratio from 65% to 55%, and employing a fallback for acyclic compounds)
to better capture modified sesterterpenoid structures.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    
    Improved Heuristics:
      • Compute the Bemis–Murcko scaffold.
      • Count carbon atoms in the scaffold and in the full molecule.
      • Count fully aromatic rings.
      • Compute the fraction of sp3-hybridized carbons and the aromatic carbon fraction.
      • Compute the ratio of scaffold carbons to full-molecule carbons.
      
    We now allow:
      - For molecules that yield a nontrivial scaffold (i.e. containing rings), a scaffold carbon count between 14 and 32 is acceptable.
      - A scaffold fraction cutoff lowered to 55% rather than 65% (to capture compounds like Zizanin-B).
      - An sp3 fraction threshold of ≥0.6 for cyclic cases.
      - For acyclic or empty scaffolds, relying on full molecule properties is preferred—with a slightly relaxed sp3 cutoff (≥0.4) and a check that the total carbon count is near a multiple of five.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a sesterterpenoid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Attempt to compute the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating scaffold: {e}"
    
    # Count carbons in the full molecule.
    full_atoms = mol.GetAtoms()
    full_mol_carbons = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6)
    
    # Count sp3 carbons in the full molecule.
    sp3_count = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6 and 
                     atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    sp3_fraction = sp3_count / full_mol_carbons if full_mol_carbons > 0 else 0
    
    # Count aromatic carbons in the full molecule.
    aromatic_carbon_count = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    aromatic_fraction = aromatic_carbon_count / full_mol_carbons if full_mol_carbons > 0 else 0
    
    # Try to compute important scaffold properties.
    # If the scaffold is essentially empty (e.g. acyclic terpenoids), we will use full molecule statistics.
    scaffold_atoms = scaffold.GetAtoms() if scaffold is not None else []
    if len(scaffold_atoms) < 3:
        scaffold_carbons = 0
        scaffold_ratio = 0
    else:
        scaffold_carbons = sum(1 for atom in scaffold_atoms if atom.GetAtomicNum() == 6)
        scaffold_ratio = scaffold_carbons / full_mol_carbons if full_mol_carbons > 0 else 0

    # Count fully aromatic rings in the scaffold.
    aromatic_ring_count = 0
    if len(scaffold_atoms) >= 3:
        ring_info = scaffold.GetRingInfo()
        for ring in ring_info.AtomRings():
            if all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_ring_count += 1

    # Check if the molecule is cyclic (i.e. has a nontrivial scaffold).
    is_cyclic = (len(scaffold_atoms) >= 3)

    # Helper: check if full molecule carbon count is near a multiple of 5 (within ±1).
    def near_mult_of_five(n):
        r = n % 5
        return r <= 1 or r >= 4

    # Primary heuristic for cyclic molecules (nontrivial Murcko scaffolds).
    if is_cyclic:
        # Allow a slightly lower scaffold carbon count (≥14 instead of 18) and lower ratio (≥55%)
        if (14 <= scaffold_carbons <= 32 and aromatic_ring_count <= 1 
            and sp3_fraction >= 0.6 and scaffold_ratio >= 0.55 and aromatic_fraction < 0.3):
            return True, (f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s), "
                          f"sp3 fraction of {sp3_fraction:.0%}, and represents {scaffold_ratio:.0%} of "
                          f"{full_mol_carbons} total carbons. These features are consistent with a sesterterpenoid backbone.")
    
    # Fallback for acyclic molecules or those with empty scaffolds 
    # (e.g. (2Z,6E,10E,14E)-geranylfarnesol might yield a nearly empty scaffold).
    if not is_cyclic:
        if (20 <= full_mol_carbons <= 60 and aromatic_fraction < 0.3 and sp3_fraction >= 0.4 
            and near_mult_of_five(full_mol_carbons)):
            return True, (f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic content "
                          f"and sp3 fraction of {sp3_fraction:.0%}. Carbon count being near a multiple of 5 "
                          f"suggests a terpenoid origin.")
    
    # As a second fallback for cyclic molecules that did not pass the stricter scaffold conditions,
    # accept if overall full-molecule statistics fall into a reasonable window and the scaffold still makes
    # up at least 55% of the carbons.
    if (22 <= full_mol_carbons <= 60 and aromatic_fraction < 0.3 and sp3_fraction >= 0.5 
        and scaffold_ratio >= 0.55 and near_mult_of_five(full_mol_carbons)):
        return True, (f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic content "
                      f"and sp3 fraction of {sp3_fraction:.0%}; scaffold ({scaffold_carbons} carbons) constitutes "
                      f"{scaffold_ratio:.0%} of the total. These values are consistent with a modified sesterterpenoid skeleton.")
    
    # If none of the conditions are met, the molecule is not considered a sesterterpenoid.
    return False, (f"Scaffold has {scaffold_carbons} carbons and represents {scaffold_ratio:.0%} of "
                   f"{full_mol_carbons} total carbons; full molecule has {aromatic_fraction:.0%} aromatic carbon "
                   f"content and sp3 fraction of {sp3_fraction:.0%}. These characteristics are not consistent with "
                   f"a typical sesterterpenoid backbone.")

# Example usage:
if __name__ == "__main__":
    # Test the function with one of the provided sesterterpenoid examples: Zizanin-B.
    test_smiles = "O=C1[C@@H]2[C@@H](C=O)CC[C@H]3[C@@](O)([C@H](CCCC(C)C)C)CC[C@@]3(C[C@@H]2[C@@](C1)(O)C)C"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)