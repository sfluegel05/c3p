"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin – any phospholipid in which the amino group of the sphingoid base 
is in amide linkage with a fatty acid and the terminal hydroxyl of that base is esterified 
to a phosphorylcholine headgroup. In addition the molecule should not contain a glycerol backbone.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    The following checks are applied:
      1. Presence of a phosphorylcholine headgroup.
         (We search for two variants of its SMARTS.)
      2. Absence of a glycerol backbone (to avoid misclassifying glycerophospholipids).
         (We search for the common glycerol fragment "OCC(O)CO".)
      3. Presence of an amide bond (the fatty acid chain is attached via amide).
      4. Presence of a sphingoid base fragment.
         (Here we require that the nitrogen attached via the amide is connected 
         to a chiral carbon bearing an –OH that is esterified to a phosphorylcholine-like fragment.)
      5. Presence of long aliphatic (fatty acid) chains.
         (At least two unrelated substructures of eight or more consecutive carbons are expected.)
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as sphingomyelin, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # 1. Check for phosphorylcholine headgroup in one of two variants.
    phos_smarts_list = [
        "P(=O)([O-])OCC[N+](C)(C)C",  # with explicit negative charge on oxygen
        "P(=O)(O)OCC[N+](C)(C)C"       # variant without explicit charge
    ]
    head_found = False
    for smarts in phos_smarts_list:
        phos_pattern = Chem.MolFromSmarts(smarts)
        if phos_pattern and mol.HasSubstructMatch(phos_pattern):
            head_found = True
            break
    if not head_found:
        return False, "Missing phosphorylcholine headgroup"
    
    # 2. Exclude glycerol backbone (common in glycerophospholipids) 
    #    using a SMARTS pattern for a glycerol fragment.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_pattern and mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone detected; likely a glycerophospholipid"
    
    # 3. Check for an amide linkage (i.e. a carbonyl attached to N).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage for fatty acid attachment"
    
    # 4. Check for a sphingoid base fragment.
    #    We require that after the amide the nitrogen is attached to a chiral carbon bearing an OH
    #    that is in turn connected to a phosphorylcholine-like fragment.
    sphingoid_smarts_patterns = [
        "[N;!R][C@@H](O)COP",   # one stereochemical assignment
        "[N;!R][C@H](O)COP"     # the other possible assignment
    ]
    sph_found = False
    for pattern_str in sphingoid_smarts_patterns:
        sph_pattern = Chem.MolFromSmarts(pattern_str)
        if sph_pattern and mol.HasSubstructMatch(sph_pattern):
            sph_found = True
            break
    if not sph_found:
        return False, "Missing sphingoid base fragment (expected N–chiral C(O)COP linkage)"
    
    # 5. Check for long aliphatic (fatty acid / sphingoid) chains.
    #    Look for at least two stretches of 8 consecutive carbon atoms.
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Missing long aliphatic chains (at least two long chain fragments expected)"
    
    return True, "Contains phosphorylcholine headgroup, amide-linked fatty acid, and a sphingoid base fragment without glycerol backbone"

# Example usage:
if __name__ == '__main__':
    test_smiles = "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
    result, reason = is_sphingomyelin(test_smiles)
    print("Result:", result)
    print("Reason:", reason)