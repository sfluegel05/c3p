"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid

A bisbenzylisoquinoline alkaloid is defined as a benzylisoquinoline alkaloid whose structure is 
built of two benzylisoquinoline units linked by ether bridges. Often, further bridging via direct 
carbon–carbon bonds or methylenedioxy groups is observed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    The approach is as follows:
      - Parse the SMILES string into an RDKit molecule.
      - Look for at least two isoquinoline-like fragments. For our purposes, we use an approximate
        SMARTS pattern for an isoquinoline core: "c1ccc2ncccc2c1".
      - Check for a bridging moiety. Typically, the two units are connected by an aromatic 
        ether bridge. In addition, methylenedioxy bridging patterns are common. We therefore 
        check for either an aromatic ether bond ([a]O[a]) or an aromatic methylenedioxy group ([a]OCO[a]).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a bisbenzylisoquinoline alkaloid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for an isoquinoline-like fragment.
    # Note: This is an approximation that looks for an aromatic system containing a nitrogen.
    isoquinoline_smarts = Chem.MolFromSmarts("c1ccc2ncccc2c1")
    iso_matches = mol.GetSubstructMatches(isoquinoline_smarts)
    
    if len(iso_matches) < 2:
        return False, f"Found {len(iso_matches)} isoquinoline-like unit(s); at least 2 required"

    # Define a SMARTS pattern for an aromatic ether bridge (an oxygen between two aromatic atoms)
    ether_bridge = Chem.MolFromSmarts("[a]O[a]")
    # Also define a pattern for a methylenedioxy group bridging two aromatic rings,
    # which is common in these alkaloids.
    md_bridge = Chem.MolFromSmarts("[a]OCO[a]")
    
    has_bridge = mol.HasSubstructMatch(ether_bridge) or mol.HasSubstructMatch(md_bridge)
    if not has_bridge:
        return False, "No bridging pattern (aromatic ether or methylenedioxy) found linking isoquinoline units"
    
    # (Optional additional filters: for example, one might check that the molecular weight or size is above a threshold.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low to be a bisbenzylisoquinoline alkaloid"
    
    return True, "Contains at least two isoquinoline-like units linked by a bridging (ether or methylenedioxy) group"

# Example test cases (the following lines can be uncommented for debugging purposes):
# test_smiles = "COc1cc2CCN(C)[C@H]3Cc4ccc(Oc5cc(C[C@H]6N(C)CCc7cc(OC)c(Oc(c1O)c23)cc67)ccc5O)cc4"  # Thalidasine example
# print(is_bisbenzylisoquinoline_alkaloid(test_smiles))