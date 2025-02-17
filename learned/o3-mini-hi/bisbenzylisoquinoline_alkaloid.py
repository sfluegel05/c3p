"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid

A bisbenzylisoquinoline alkaloid is defined as a benzylisoquinoline alkaloid whose structures are 
built up of two benzylisoquinoline units linked by ether bridges. Various additional bridging patterns 
(direct carbon–carbon bonds or methylenedioxy groups) can be observed.

Our revised criteria are:
  - The molecule must be a valid SMILES and have molecular weight >= 500 Da.
  - It must contain at least two distinct isoquinoline-like substructures. We now define two looser 
    SMARTS patterns to capture:
      • an aromatic isoquinoline: "c1ccc2c(c1)cccn2"
      • a tetrahydroisoquinoline-like unit: "c1ccc2CCNC2c1"
    Unique matches from either pattern are combined.
  - At least one bridging motif connecting aromatic portions must be present, either:
      • an aromatic ether bridge ([a]O[a]) or
      • a methylenedioxy bridge ([a]OCO[a]).
  - Finally, as benzylisoquinoline units normally carry a nitrogen, the molecule should contain at least 2 nitrogen atoms.
  
If all criteria are met, the function returns (True, <reason>), otherwise (False, <reason>).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    
    The function checks:
      - Valid SMILES string and molecular weight >= 500 Da.
      - Presence of at least two unique isoquinoline-like substructures. Two patterns are used:
          (a) Aromatic isoquinoline, SMARTS: "c1ccc2c(c1)cccn2"
          (b) Tetrahydroisoquinoline-like fragment, SMARTS: "c1ccc2CCNC2c1"
      - Presence of a bridging motif between aromatic portions: either an aromatic ether bridge ([a]O[a]) 
        or a methylenedioxy bridge ([a]OCO[a]).
      - At least 2 nitrogen atoms in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule qualifies as a bisbenzylisoquinoline alkaloid, 
                     and False otherwise; the second element is an explanation.
    """
    # Convert SMILES to molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (must be >= 500 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a bisbenzylisoquinoline alkaloid"
    
    # Ensure there are at least 2 nitrogen atoms in the molecule.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, f"Not enough nitrogen atoms ({n_count}); at least 2 are expected in bisbenzylisoquinoline alkaloids"
    
    # Define updated SMARTS patterns for isoquinoline-like substructures.
    # Pattern (a): Aromatic isoquinoline (benzene fused with a pyridine-like ring).
    arom_iso_smarts = "c1ccc2c(c1)cccn2"
    arom_iso = Chem.MolFromSmarts(arom_iso_smarts)
    if arom_iso is None:
        return False, "Error creating aromatic isoquinoline SMARTS pattern"
    
    # Pattern (b): Tetrahydroisoquinoline-like unit (benzene fused with a saturated N–containing ring).
    tetra_iso_smarts = "c1ccc2CCNC2c1"
    tetra_iso = Chem.MolFromSmarts(tetra_iso_smarts)
    if tetra_iso is None:
        return False, "Error creating tetrahydroisoquinoline SMARTS pattern"
    
    # Find all substructure matches for both patterns.
    arom_matches = mol.GetSubstructMatches(arom_iso, uniquify=True)
    tetra_matches = mol.GetSubstructMatches(tetra_iso, uniquify=True)
    
    # Combine matches uniquely based on the set of atom indices.
    unique_matches = set()
    for match in arom_matches:
        unique_matches.add(frozenset(match))
    for match in tetra_matches:
        unique_matches.add(frozenset(match))
        
    if len(unique_matches) < 2:
        return False, f"Found only {len(unique_matches)} isoquinoline-like substructure(s); at least 2 are required"

    # Look for a bridging motif connecting aromatic fragments.
    # Pattern for aromatic ether bridge: oxygen connected to two aromatic atoms.
    ether_bridge_smarts = "[a]O[a]"
    ether_bridge = Chem.MolFromSmarts(ether_bridge_smarts)
    if ether_bridge is None:
        return False, "Error creating aromatic ether bridge SMARTS pattern"
    
    # Pattern for methylenedioxy bridge: -OCO- bridge between two aromatic atoms.
    md_bridge_smarts = "[a]OCO[a]"
    md_bridge = Chem.MolFromSmarts(md_bridge_smarts)
    if md_bridge is None:
        return False, "Error creating methylenedioxy bridge SMARTS pattern"
    
    has_bridge = mol.HasSubstructMatch(ether_bridge) or mol.HasSubstructMatch(md_bridge)
    if not has_bridge:
        return False, "No bridging pattern (aromatic ether or methylenedioxy) found linking aromatic portions"
    
    return True, ("Molecule has molecular weight {:.1f} Da, contains at least 2 benzylisoquinoline-like substructures, "
                  "has at least 2 nitrogen atoms, and exhibits a bridging pattern consistent with two benzylisoquinoline units"
                  .format(mol_wt))
    
# Example usage (for testing/debugging):
# test_smiles = "COc1ccc2C[C@@H]3N(C)[C@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(O)c(Oc1c2)cc7)c67)cc3)c2cc1OC"  # Example: Thalidasine
# print(is_bisbenzylisoquinoline_alkaloid(test_smiles))