"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin – any phospholipid in which the amino group of a sphingoid base 
is in amide linkage with a fatty acid and the terminal hydroxyl of that base is esterified 
to a phosphorylcholine headgroup. In addition the molecule should not contain a glycerol backbone.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    The following checks are applied:
      1. Verify the presence of a phosphorylcholine headgroup.
         (Two variants are considered to account for explicit vs implicit charges.)
      2. Exclude molecules that contain a glycerol backbone (common in glycerophospholipids).
      3. Check for the presence of an amide linkage.
      4. Look for a sphingoid base fragment.
         Instead of requiring an explicit hydroxyl on the chiral carbon,
         we search for an amide-linked nitrogen attached to a chiral carbon that is directly 
         bonded to a "COP" fragment (capturing the headgroup esterification).
      5. Ensure the presence of long aliphatic (fatty acid/sphingoid) chains.
         At least two substructures of eight or more consecutive carbons are expected.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as sphingomyelin, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # 1. Check for phosphorylcholine headgroup (two variants)
    pc_smarts_list = [
        "P(=O)([O-])OCC[N+](C)(C)C",  # with explicit negative charge
        "P(=O)(O)OCC[N+](C)(C)C"       # variant without explicit charge
    ]
    pc_found = False
    for smarts in pc_smarts_list:
        pc_pattern = Chem.MolFromSmarts(smarts)
        if pc_pattern and mol.HasSubstructMatch(pc_pattern):
            pc_found = True
            break
    if not pc_found:
        return False, "Missing phosphorylcholine headgroup"
    
    # 2. Exclude glycerol backbone (common in glycerophospholipids).
    #    Using a SMARTS pattern for the glycerol fragment "OCC(O)CO"
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_pattern and mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone detected; likely a glycerophospholipid"
    
    # 3. Check for an amide linkage (C(=O)N pattern).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage for fatty acid attachment"
    
    # 4. Check for a sphingoid base fragment.
    #    Instead of requiring an explicit (O) on the chiral carbon,
    #    we look for an N (not in a ring) connected to a chiral carbon that is bonded to a COP fragment.
    sphingoid_patterns = [
        "[N;!R][C@@H]COP",   # chiral carbon with phosphorylcholine fragment (one stereochemistry)
        "[N;!R][C@H]COP"     # the other possible stereochemical assignment
    ]
    sph_found = False
    for pat_str in sphingoid_patterns:
        sph_pattern = Chem.MolFromSmarts(pat_str)
        if sph_pattern and mol.HasSubstructMatch(sph_pattern):
            sph_found = True
            break
    if not sph_found:
        return False, "Missing sphingoid base fragment (expected N–chiral C(COP) linkage)"
    
    # 5. Check for long aliphatic chains: at least two stretches of 8 consecutive carbons.
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(chain_matches) < 2:
        return False, "Missing long aliphatic chains (require at least two stretches of 8+ carbons)"
    
    return True, "Contains phosphorylcholine headgroup, amide-linked fatty acid, and a sphingoid base fragment without glycerol backbone"

# Example usage:
if __name__ == '__main__':
    test_smiles = "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
    result, reason = is_sphingomyelin(test_smiles)
    print("Result:", result)
    print("Reason:", reason)