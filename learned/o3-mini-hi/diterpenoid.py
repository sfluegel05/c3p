"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: Diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 skeleton 
of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms 
(generally methyl groups).

Heuristic criteria used in this function:
  - Allowed elements are C, H and O only (this removes many false‐positives that include e.g. nitrogen, halogens, etc.).
  - The number of carbon atoms is expected to be roughly between 15 and 33.
  - The molecular weight is expected to be in the range 220–800 Da.
  - Many diterpenoids feature at least one non‐aromatic (aliphatic) ring.
If these criteria are met, the molecule is considered likely to be a diterpenoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string using heuristic checks.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a diterpenoid, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Criterion 1: Allowed Elements. Diterpenoids are biosynthesized from isoprene units and typically contain only C, H, and O.
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()} which is not typical for diterpenoids."

    # Criterion 2: Carbon count.
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if nC < 15:
        return False, f"Too few carbon atoms ({nC}) for a diterpenoid skeleton."
    if nC > 33:
        return False, f"Too many carbon atoms ({nC}); exceeds expected diterpenoid range."

    # Criterion 3: Molecular weight check.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 220 or mw > 800:
        return False, f"Molecular weight {mw:.1f} Da is outside the typical range (220-800 Da) for diterpenoids."

    # Criterion 4: Ring system analysis – many diterpenoids contain at least one non‐aromatic (aliphatic) ring.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 1:
        return False, "No cyclic ring system detected; many diterpenoids are cyclic."
    # Count number of rings that are not entirely aromatic.
    non_aromatic_rings = 0
    for ring in ring_info.AtomRings():
        # Check every atom in the ring to see if it is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            non_aromatic_rings += 1
    if non_aromatic_rings < 1:
        return False, "All ring systems are aromatic; diterpenoids usually have aliphatic ring(s)."

    # If all heuristic checks are passed, we classify the molecule as likely a diterpenoid.
    return True, (f"Likely diterpenoid: {nC} carbons, {mw:.1f} Da molecular weight, " +
                  f"{n_rings} ring(s) detected with at least {non_aromatic_rings} non-aromatic ring(s).")

# Example use (these lines can be removed or commented out in production):
if __name__ == "__main__":
    test_smiles = "C[C@H](CCC1OC1(C)C)[C@@H]1CCC(C)=C2CC=C(C)[C@H]2[C@@H]1O"  # Acutilol A
    result, reason = is_diterpenoid(test_smiles)
    print(result, reason)