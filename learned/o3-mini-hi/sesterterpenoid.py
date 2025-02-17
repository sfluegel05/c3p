"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (a C25 precursor) including compounds 
where the C25 skeleton has been rearranged or partly truncated.
This improved script uses a combination of Bemis–Murcko scaffold analysis and full‐molecule 
statistics (carbon count, aromaticity, sp3 fraction, and scaffold-to-full carbon ratio) while 
relaxing some thresholds in order to capture more modified sesterterpenoid structures and reduce 
false negatives. At the same time, a fallback for mostly acyclic cases is used with a narrow full-molecule 
carbon range to help minimize false positives.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    
    Improved Heuristics:
      • Compute the Bemis–Murcko scaffold.
      • Count carbon atoms in the scaffold and in the full molecule.
      • Evaluate the number of fully aromatic rings.
      • Compute the sp3 carbon fraction and aromatic carbon fraction for the full molecule.
      • Compute the ratio of scaffold carbons to full-molecule carbons.
      
    We now allow:
      - For molecules with a nontrivial scaffold (i.e. cyclic), a scaffold carbon count threshold of 14–32,
        a minimum scaffold-to-full ratio of 50% (instead of 55%), and sp3 fraction of ≥0.6.
      - For acyclic molecules (or those with an empty scaffold) a fallback based solely on full-molecule properties:
        full carbon count must lie between 22 and 35 (this helps exclude simple fatty acids) plus a relaxed
        sp3 fraction (≥0.4) and an overall carbon count near a multiple of 5.
      
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
    
    # Attempt to compute the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating scaffold: {e}"
    
    # Count carbons in the full molecule.
    full_atoms = mol.GetAtoms()
    full_mol_carbons = sum(1 for atom in full_atoms if atom.GetAtomicNum() == 6)
    if full_mol_carbons == 0:
        return False, "No carbon atoms found in molecule."
    
    # Count sp3 carbons in full molecule.
    sp3_count = sum(1 for atom in full_atoms if (atom.GetAtomicNum() == 6 and 
                                                   atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3))
    sp3_fraction = sp3_count / full_mol_carbons
    
    # Count aromatic carbons in full molecule.
    aromatic_carbons = sum(1 for atom in full_atoms if (atom.GetAtomicNum() == 6 and atom.GetIsAromatic()))
    aromatic_fraction = aromatic_carbons / full_mol_carbons
    
    # Compute scaffold statistics.
    scaffold_atoms = scaffold.GetAtoms() if scaffold is not None else []
    # If the scaffold is too small (for example, acyclic molecules may yield 0 or 1 atoms),
    # treat it as "empty" by our criteria.
    if len(scaffold_atoms) < 3:
        scaffold_carbons = 0
        scaffold_ratio = 0
    else:
        scaffold_carbons = sum(1 for atom in scaffold_atoms if atom.GetAtomicNum() == 6)
        scaffold_ratio = scaffold_carbons / full_mol_carbons

    # Count fully aromatic rings in the scaffold.
    aromatic_ring_count = 0
    if len(scaffold_atoms) >= 3:
        ring_info = scaffold.GetRingInfo()
        for ring in ring_info.AtomRings():
            if all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_ring_count += 1

    # Determine if the molecule is cyclic (i.e. nontrivial scaffold exists).
    is_cyclic = (len(scaffold_atoms) >= 3)

    # Helper: check if a number is "near" a multiple of 5 (within ±1).
    def near_mult_of_five(n):
        r = n % 5
        return r <= 1 or r >= 4

    # Primary heuristic for cyclic molecules.
    if is_cyclic:
        # Check if scaffold properties and overall molecular metrics fall in the desired range.
        # Relaxed thresholds:
        #   - scaffold_carbons: between 14 and 32
        #   - scaffold_ratio: at least 50%
        #   - sp3 fraction: at least 60%
        #   - aromatic rings in scaffold: at most 1
        #   - overall aromatic carbon fraction: below 30%
        if (14 <= scaffold_carbons <= 32 and aromatic_ring_count <= 1 and
            sp3_fraction >= 0.6 and scaffold_ratio >= 0.50 and aromatic_fraction < 0.3):
            return True, (f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s), "
                          f"sp3 fraction of {sp3_fraction:.0%}, and represents {scaffold_ratio:.0%} of "
                          f"{full_mol_carbons} total carbons. These features are consistent with a sesterterpenoid backbone.")
    
    # Fallback for molecules lacking a meaningful scaffold (acyclic)
    # Here we want to avoid long, unbranched fatty acids.
    if not is_cyclic:
        if (22 <= full_mol_carbons <= 35 and aromatic_fraction < 0.3 and sp3_fraction >= 0.4 and
            near_mult_of_five(full_mol_carbons)):
            return True, (f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic content "
                          f"and sp3 fraction of {sp3_fraction:.0%}. Carbon count being near a multiple of 5 suggests "
                          f"a terpenoid origin.")
    
    # Second fallback for cyclic molecules that did not meet the stricter scaffold criteria:
    # Use full-molecule metrics with a slightly relaxed sp3 cutoff and scaffold ratio lowered to 50%.
    if (22 <= full_mol_carbons <= 60 and aromatic_fraction < 0.3 and sp3_fraction >= 0.5 and
        scaffold_ratio >= 0.50 and near_mult_of_five(full_mol_carbons)):
        return True, (f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic content "
                      f"and sp3 fraction of {sp3_fraction:.0%}; scaffold ({scaffold_carbons} carbons) constitutes "
                      f"{scaffold_ratio:.0%} of the total. These values are consistent with a modified sesterterpenoid skeleton.")
    
    # If none of these conditions are met, classify as not a sesterterpenoid.
    return False, (f"Scaffold has {scaffold_carbons} carbons (making up {scaffold_ratio:.0%} of {full_mol_carbons} total carbons); "
                   f"full molecule has {aromatic_fraction:.0%} aromatic carbons and sp3 fraction of {sp3_fraction:.0%}. "
                   f"These characteristics are not consistent with a typical sesterterpenoid backbone.")

# Example usage:
if __name__ == "__main__":
    # Test the function with one of the reference sesterterpenoid examples: Zizanin-B.
    test_smiles = "O=C1[C@@H]2[C@@H](C=O)CC[C@H]3[C@@](O)([C@H](CCCC(C)C)C)CC[C@@]3(C[C@@H]2[C@@](C1)(O)C)C"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)