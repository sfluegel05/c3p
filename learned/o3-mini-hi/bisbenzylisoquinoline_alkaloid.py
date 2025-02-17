"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid

A bisbenzylisoquinoline alkaloid is defined as a benzylisoquinoline alkaloid whose structures are 
built from two benzylisoquinoline units linked by ether bridges. Often, further bridging by direct 
carbon–carbon bonds or methylenedioxy groups is observed.

This script uses the following criteria:
  - Valid SMILES string that generates an RDKit molecule.
  - Molecular weight >= 500 Da.
  - Presence of at least 2 benzylisoquinoline-like substructures. This is determined by 
    combining matches from an aromatic isoquinoline SMARTS ("c1ccc2nc(ccc2c1)") and 
    a tetrahydroisoquinoline SMARTS ("c1ccc2CCN(C2)c1").
  - A bridging motif between aromatic portions, meaning an aromatic ether bridge ([a]O[a]) 
    or a methylenedioxy bridge ([a]OCO[a]).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    The approach is as follows:
      - Convert the SMILES to an RDKit molecule.
      - Check the molecular weight (should be >= 500 Da).
      - Look for at least two distinct benzylisoquinoline-like substructures. This is done by detecting 
        either aromatic isoquinoline (SMARTS: "c1ccc2nc(ccc2c1)") or tetrahydroisoquinoline-like units
        (SMARTS: "c1ccc2CCN(C2)c1").
      - Verify that a bridging motif exists linking aromatic parts: either an aromatic ether bridge ([a]O[a])
        or a methylenedioxy pattern ([a]OCO[a]).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element is True if the molecule is classified as a bisbenzylisoquinoline alkaloid,
                     and False otherwise; the second element is a reason string.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a bisbenzylisoquinoline alkaloid"
    
    # Define SMARTS patterns for isoquinoline-like substructures.
    # Pattern 1: Fully aromatic isoquinoline: benzene fused to pyridine ring.
    isoquinoline_smarts = "c1ccc2nc(ccc2c1)"
    arom_iso = Chem.MolFromSmarts(isoquinoline_smarts)
    if arom_iso is None:
        return False, "Error creating aromatic isoquinoline SMARTS pattern"
    
    # Pattern 2: Tetrahydroisoquinoline-like unit: benzene fused to a saturated N-containing ring.
    # This pattern allows for a methlylated nitrogen as well.
    tetra_iso_smarts = "c1ccc2CCN(C2)c1"
    tetra_iso = Chem.MolFromSmarts(tetra_iso_smarts)
    if tetra_iso is None:
        return False, "Error creating tetrahydroisoquinoline SMARTS pattern"
    
    # Find all substructure matches from both patterns.
    arom_matches = mol.GetSubstructMatches(arom_iso, uniquify=True)
    tetra_matches = mol.GetSubstructMatches(tetra_iso, uniquify=True)
    
    # Combine matches uniquely (based on the set of atom indices for each match).
    unique_matches = set()
    for match in arom_matches:
        unique_matches.add(frozenset(match))
    for match in tetra_matches:
        unique_matches.add(frozenset(match))
    
    if len(unique_matches) < 2:
        return False, f"Found only {len(unique_matches)} isoquinoline-like substructure(s); at least 2 are required"
    
    # Look for a bridging motif linking the aromatic units.
    # Aromatic ether bridge: oxygen connected to two aromatic atoms.
    ether_bridge_smarts = "[a]O[a]"
    ether_bridge = Chem.MolFromSmarts(ether_bridge_smarts)
    if ether_bridge is None:
        return False, "Error creating aromatic ether bridge SMARTS pattern"
    
    # Methylenedioxy bridge: pattern where oxygen is in a -OCO- group bridging two aromatic atoms.
    md_bridge_smarts = "[a]OCO[a]"
    md_bridge = Chem.MolFromSmarts(md_bridge_smarts)
    if md_bridge is None:
        return False, "Error creating methylenedioxy bridge SMARTS pattern"
    
    has_bridge = mol.HasSubstructMatch(ether_bridge) or mol.HasSubstructMatch(md_bridge)
    if not has_bridge:
        return False, "No bridging pattern (aromatic ether or methylenedioxy) found linking aromatic portions"
    
    return True, ("Molecule has molecular weight {:.1f} Da, contains at least 2 benzylisoquinoline-like substructures, and "
                  "exhibits a bridging pattern consistent with two benzylisoquinoline units".format(mol_wt))

# Example usage for debugging:
# test_smiles = "COc1ccc2C[C@@H]3N(C)[C@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(OC)c(Oc1c2)cc7)c67)cc3)c2cc1OC"  # Example: Thalidasine
# print(is_bisbenzylisoquinoline_alkaloid(test_smiles))