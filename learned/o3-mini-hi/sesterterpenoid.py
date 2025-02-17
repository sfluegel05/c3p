"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (a C25 precursor), 
including compounds where the C25 skeleton has been rearranged or partly truncated.
This script uses multiple heuristics based on the Bemis–Murcko scaffold (its carbon count,
aromatic ring count, and sp3 character) as well as full-molecule statistics to decide on the class.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    
    Heuristics:
      1. Compute the Bemis–Murcko scaffold, count its carbon atoms,
         count fully aromatic rings in the scaffold, and assess its sp3 character.
         We expect a terpenoid core (from a C25 precursor) to have roughly 20–30 carbons,
         to be largely aliphatic (≤1 aromatic ring), and to have a high fraction of sp3 carbons.
      2. Failing that, we look at the full molecule:
         we require an overall carbon count roughly in the range (22–40) with low aromatic fraction (<30%)
         and a moderate sp3 fraction.
      
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
    
    # Compute the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating scaffold: {e}"
    if scaffold is None:
        return False, "Could not compute Murcko scaffold for molecule."
    
    # Count carbon atoms in the scaffold.
    scaffold_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count fully aromatic rings in the scaffold.
    aromatic_ring_count = 0
    ri = scaffold.GetRingInfo()
    for ring in ri.AtomRings():
        if all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    # Count sp3-hybridized carbons in full molecule.
    full_mol_atoms = mol.GetAtoms()
    full_mol_carbons = sum(1 for atom in full_mol_atoms if atom.GetAtomicNum() == 6)
    sp3_count = sum(1 for atom in full_mol_atoms if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    sp3_fraction = (sp3_count / full_mol_carbons) if full_mol_carbons > 0 else 0

    # Count aromatic carbons in the full molecule.
    aromatic_carbon_count = sum(1 for atom in full_mol_atoms if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    aromatic_fraction = (aromatic_carbon_count / full_mol_carbons) if full_mol_carbons > 0 else 0

    # Primary heuristic: Evaluate the scaffold.
    # We expect a sesterterpenoid core (modified C25) to have roughly 20-30 carbons,
    # minimal aromatic rings (≤1), and high sp3 fraction (≥0.6).
    if 20 <= scaffold_carbons <= 30 and aromatic_ring_count <= 1 and sp3_fraction >= 0.6:
        return True, (f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s) "
                      f"and sp3 fraction of {sp3_fraction:.0%}, consistent with a sesterterpenoid backbone.")
    
    # Fallback heuristic: Evaluate the full molecule.
    # Allow a slightly wider range for full molecule carbon count (22-40) and require low aromatic fraction (<30%)
    # and decent sp3 character (≥0.5).
    if 22 <= full_mol_carbons <= 40 and aromatic_fraction < 0.3 and sp3_fraction >= 0.5:
        return True, (f"Full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic carbon content "
                      f"and sp3 fraction of {sp3_fraction:.0%}, consistent with a modified sesterterpenoid skeleton.")
    
    # Otherwise, do not classify as a sesterterpenoid.
    return False, (f"Scaffold has {scaffold_carbons} carbons with {aromatic_ring_count} aromatic ring(s); "
                   f"full molecule has {full_mol_carbons} carbons with {aromatic_fraction:.0%} aromatic content "
                   f"and sp3 fraction of {sp3_fraction:.0%}, which is not consistent with a typical sesterterpenoid.")

# Example usage:
if __name__ == "__main__":
    # Using (2Z,6E,10E,14E)-geranylfarnesol as an example (a known sesterterpenoid precursor).
    test_smiles = "C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)/C(=C\\CO)/C"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)