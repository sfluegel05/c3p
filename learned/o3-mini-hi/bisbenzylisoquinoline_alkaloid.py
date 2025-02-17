"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid

A bisbenzylisoquinoline alkaloid is defined as a benzylisoquinoline alkaloid whose structure is 
built from two benzylisoquinoline units linked by ether bridges. Often, further bridging by direct 
carbon–carbon bonds or methylenedioxy groups is observed.

The improved criteria used in this script are:
  - Valid SMILES that generate an RDKit molecule.
  - A molecular weight above 500 Da.
  - The presence of at least two isoquinoline-like substructures.
  - A bridging moiety between aromatic parts as seen by an aromatic ether bridge ([a]O[a]) or a 
    methylenedioxy bridging pattern ([a]OCO[a]).
Note: This heuristic approach uses SMARTS patterns as a proxy for a true substructure classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    The approach is as follows:
      - Convert the SMILES to an RDKit molecule.
      - Check the molecular weight (should be >= 500 Da for these relatively heavy dimers).
      - Look for at least two distinct isoquinoline substructures (representing the benzylisoquinoline units).
      - Verify that a bridging motif exists linking aromatic portions (either an aromatic ether bridge or
        a methylenedioxy pattern).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element is True if classified as a bisbenzylisoquinoline alkaloid,
                     and False otherwise; the second element is a reason string.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a bisbenzylisoquinoline alkaloid"
    
    # Define a SMARTS pattern to detect an isoquinoline-like ring.
    # This pattern looks for a fused benzene and pyridine ring: e.g., c1ccc2nc(ccc2c1)
    isoquinoline_smarts = "c1ccc2nc(ccc2c1)"
    isoquinoline_pattern = Chem.MolFromSmarts(isoquinoline_smarts)
    if isoquinoline_pattern is None:
        return False, "Error creating isoquinoline SMARTS pattern"
    
    iso_matches = mol.GetSubstructMatches(isoquinoline_pattern, uniquify=True)
    # We need at least 2 distinct isoquinoline substructures.
    if len(iso_matches) < 2:
        return False, f"Found only {len(iso_matches)} isoquinoline substructure(s); at least 2 are required"
    
    # Look for a bridging motif.
    # Aromatic ether bridge: an oxygen connected to two aromatic atoms.
    ether_bridge_smarts = "[a]O[a]"
    ether_bridge = Chem.MolFromSmarts(ether_bridge_smarts)
    # Methylenedioxy bridge: pattern [a]OCO[a]
    md_bridge_smarts = "[a]OCO[a]"
    md_bridge = Chem.MolFromSmarts(md_bridge_smarts)
    
    has_bridge = mol.HasSubstructMatch(ether_bridge) or mol.HasSubstructMatch(md_bridge)
    if not has_bridge:
        return False, "No bridging pattern (aromatic ether or methylenedioxy) found linking aromatic portions"
    
    return True, ("Molecule has molecular weight {:.1f} Da, contains at least 2 isoquinoline-like substructures, and "
                  "exhibits a bridging pattern consistent with two benzylisoquinoline units".format(mol_wt))

# Example usage (for debugging purposes, uncomment the following lines):
# test_smiles = "COc1ccc2C[C@@H]3N(C)[C@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(OC)c(Oc1c2)cc7)c67)cc3)c2cc1OC"  # Example: Thalidasine
# print(is_bisbenzylisoquinoline_alkaloid(test_smiles))