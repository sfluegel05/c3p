"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: Diterpenoid
Definition: Any terpenoid derived from a diterpene. The term includes compounds in which the C20 skeleton 
of the parent diterpene has been rearranged or modified by the removal of one or more skeletal atoms 
(generally methyl groups).

Heuristic criteria used in this function:
  - The molecule should have a number of carbon atoms roughly in the range of 15 to 25.
    (The parent diterpene has 20 carbons, but modifications may remove or occasionally add carbons.)
  - The molecular weight is expected to fall in approximately 250 to 600 Da.
  - Many diterpenoids feature cyclic ring systems.
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
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count carbon atoms (atomic number 6)
    nC = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count heavy atoms (non-hydrogen)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    nHeavy = len(heavy_atoms)
    
    # Heuristic 1: Check if the number of carbons is in a plausible range.
    # Diterpenes have a C20 skeleton, but modifications may cause small deviations.
    if nC < 15:
        return False, f"Too few carbon atoms ({nC}) for a diterpenoid skeleton."
    if nC > 25:
        # If there are many extra carbons, check the fraction relative to all heavy atoms.
        ratio = nC / nHeavy if nHeavy > 0 else 0
        if ratio < 0.5:
            return False, f"Carbon fraction too low ({ratio:.2f}) for a diterpenoid skeleton."
    
    # Heuristic 2: Molecular weight check (approximate range for many diterpenoids)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250 or mw > 600:
        return False, f"Molecular weight {mw:.1f} Da is outside the typical range (250-600 Da) for diterpenoids."
    
    # Heuristic 3: Most diterpenoids are cyclic.
    ring_info = mol.GetRingInfo()
    nRings = ring_info.NumRings()
    if nRings < 1:
        return False, "No cyclic ring system detected; many diterpenoids are cyclic."
    
    # If all heuristic checks are passed, classify as diterpenoid.
    return True, (f"Likely diterpenoid: {nC} carbons, {mw:.1f} Da molecular weight, and {nRings} ring(s) detected.")