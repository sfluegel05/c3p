"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid

A bisbenzylisoquinoline alkaloid is defined as a benzylisoquinoline alkaloid whose structure is 
built from two benzylisoquinoline units linked by ether bridges. Often, further bridging by direct 
carbon–carbon bonds or methylenedioxy groups is observed.
The criteria in this script are:
  - Valid SMILES that generates a molecule.
  - A molecular weight above 500 Da (most true bisbenzylisoquinolines are relatively heavy dimers).
  - The presence of at least two aromatic nitrogen atoms (assuming each benzylisoquinoline unit contains a nitrogen).
  - A bridging moiety between aromatic parts as seen by an aromatic ether bridge ([a]O[a]) or a 
    methylenedioxy bridging pattern ([a]OCO[a]).
Note: This is an approximate method that uses heuristic SMARTS patterns.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    The approach is as follows:
      - Convert the SMILES to an RDKit molecule.
      - Check molecular weight (should be >= 500 Da for sizable dimeric alkaloids).
      - Check that there are at least two aromatic nitrogen atoms.
      - Look for a bridging pattern defined as an aromatic ether ([a]O[a]) or a methylenedioxy ([a]OCO[a])
        connecting aromatic substructures.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a bisbenzylisoquinoline alkaloid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a bisbenzylisoquinoline alkaloid"
    
    # Count aromatic nitrogen atoms.
    # We loop over all atoms; if the atom is nitrogen and is aromatic, we count it.
    aromatic_nitrogen_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic():
            aromatic_nitrogen_count += 1
    if aromatic_nitrogen_count < 2:
        return False, f"Found only {aromatic_nitrogen_count} aromatic nitrogen atom(s); at least 2 required"
    
    # Look for a bridging motif.
    # Define SMARTS for an aromatic ether bridge: an oxygen bonded to two aromatic atoms.
    ether_bridge_smarts = "[a]O[a]"
    ether_bridge = Chem.MolFromSmarts(ether_bridge_smarts)
    # Also, a methylenedioxy bridging pattern is common: two oxygen atoms with a central carbon bridging.
    md_bridge_smarts = "[a]OCO[a]"
    md_bridge = Chem.MolFromSmarts(md_bridge_smarts)
    
    has_bridge = mol.HasSubstructMatch(ether_bridge) or mol.HasSubstructMatch(md_bridge)
    if not has_bridge:
        return False, "No bridging pattern (aromatic ether or methylenedioxy) found linking aromatic portions"
    
    return True, ("Molecule has molecular weight {:.1f} Da, contains at least 2 aromatic nitrogen(s), and "
                  "exhibits a bridging pattern consistent with two benzylisoquinoline units".format(mol_wt))

# Example usage (for debugging purposes, uncomment the following lines):
# test_smiles = "COc1ccc2C[C@@H]3N(C)[C@H](Cc3ccc(Oc4cc(C[C@H]5N(C)CCc6cc(OC)c(Oc1c2)cc7)c67)cc3)c2cc1OC"  # Example: Thalidasine
# print(is_bisbenzylisoquinoline_alkaloid(test_smiles))