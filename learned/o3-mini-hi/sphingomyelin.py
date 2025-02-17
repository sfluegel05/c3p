"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin – any phospholipid in which the amino group of a sphingoid base 
is in amide linkage with a fatty acid while the terminal hydroxyl of that base is esterified 
to a phosphorylcholine headgroup. Also, the molecule should not contain a glycerol backbone and 
should be acyclic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    The checks applied are:
      1. Presence of a phosphorylcholine headgroup (using two SMARTS variants to capture charge differences).
      2. Exclusion of a glycerol backbone (to avoid glycerophospholipids).
      3. Presence of an amide linkage (C(=O)N) for the fatty acid attachment.
      4. Presence of a sphingoid base fragment – we look for an N (not in a ring) bound to a chiral carbon 
         bearing a hydroxyl and attached to a “COP” fragment (the linkage to the headgroup).
      5. Verification that the molecule does NOT contain any ring structures (most sphingomyelin molecules are acyclic).
      6. At least one long uninterrupted aliphatic chain (8 or more consecutive carbons) is required.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as sphingomyelin, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # 1. Check for phosphorylcholine headgroup.
    #    We use two SMARTS patterns to account for explicit negative charges or not.
    pc_smarts_list = [
        "P(=O)([O-])OCC[N+](C)(C)C",
        "P(=O)(O)OCC[N+](C)(C)C"
    ]
    pc_found = False
    for smarts in pc_smarts_list:
        pc_pattern = Chem.MolFromSmarts(smarts)
        if pc_pattern and mol.HasSubstructMatch(pc_pattern):
            pc_found = True
            break
    if not pc_found:
        return False, "Missing phosphorylcholine headgroup"
    
    # 2. Exclude glycerol backbone.
    #    If the common glycerol fragment "OCC(O)CO" is found, then this is likely a glycerophospholipid.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_pattern and mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone detected; likely a glycerophospholipid"
    
    # 3. Check for an amide linkage. The presence of C(=O)N is required.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage for fatty acid attachment"
    
    # 4. Check for a sphingoid base fragment.
    # We search for an N (not in a ring) attached to a chiral carbon that is bonded to a COP fragment.
    sphingoid_patterns = [
        "[N;!R][C@@H](O)COP",  # one stereochemistry
        "[N;!R][C@H](O)COP"     # the other possibility
    ]
    sph_found = False
    for pat_str in sphingoid_patterns:
        sph_pattern = Chem.MolFromSmarts(pat_str)
        if sph_pattern and mol.HasSubstructMatch(sph_pattern):
            sph_found = True
            break
    if not sph_found:
        return False, "Missing sphingoid base fragment (expected N–chiral C(O)COP linkage)"
    
    # 5. Exclude molecules with any ring structures.
    # Sphingomyelins are typically acyclic in their sphingoid and fatty acid parts.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures, unlikely to be sphingomyelin"
    
    # 6. Check for long aliphatic chain(s):
    # A chain of 8 or more consecutive carbons is expected.
    # (Relaxed requirement: at least one occurrence.)
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(chain_matches) < 1:
        return False, "Missing long aliphatic chain (require at least one stretch of 8+ carbons)"
    
    return True, "Contains phosphorylcholine headgroup, amide-linked fatty acid, and a sphingoid base fragment without glycerol backbone and rings; long aliphatic chain present"


# Example usage:
if __name__ == '__main__':
    # Test with one of the provided sphingomyelin SMILES as an example.
    test_smiles = "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
    result, reason = is_sphingomyelin(test_smiles)
    print("Result:", result)
    print("Reason:", reason)