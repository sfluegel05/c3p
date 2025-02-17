"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin – any phospholipid where the amino group of a sphingoid base 
is in amide linkage with one of several fatty acids and the terminal hydroxyl of the 
sphingoid base is esterified to a phosphorylcholine headgroup.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    Requirements:
      - It has a phosphorylcholine headgroup.
      - It contains an amide linkage (for fatty acid attachment).
      - It contains a sphingoid base fragment that can be identified by the pattern:
        a fatty acyl carbonyl attached to a nitrogen and then two consecutive chiral carbons,
        the second of which bears an OH.
      - It also has a long aliphatic (fatty acid) chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a sphingomyelin, False otherwise.
        str: Explanation of the classification decision.
    """
    
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Check for the phosphorylcholine headgroup with two variants.
    phos_smarts_list = [
        "P(=O)([O-])OCC[N+](C)(C)C",   # with explicit negative charge on oxygen
        "P(=O)(O)OCC[N+](C)(C)C"        # without explicit charge info
    ]
    headgroup_found = False
    for pattern in phos_smarts_list:
        phos_pattern = Chem.MolFromSmarts(pattern)
        if phos_pattern and mol.HasSubstructMatch(phos_pattern):
            headgroup_found = True
            break
    if not headgroup_found:
        return False, "Missing phosphorylcholine headgroup"
    
    # Check for amide bond linking a fatty acid through a carbonyl attached to N.
    # This will match for example the pattern "CCCC(=O)N" (the fatty acid part plus amide nitrogen).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage for fatty acid attachment"
    
    # Check for the sphingoid base.
    # Instead of simply searching for "N[C@H](O)", we now look for a two‐carbon fragment
    # with defined chirality where the first carbon follows a carbonyl (the acyl part) and
    # the second carries a hydroxyl.
    #
    # We use two patterns to allow either chiral assignment:
    sph_base_smarts_patterns = [
        "[CX3](=O)N[C@@H][C@H](O)",  # pattern 1: carbonyl, then N, then chiral carbon followed by chiral OH‐bearing carbon.
        "[CX3](=O)N[C@H][C@@H](O)"   # pattern 2: swapped chirality labels.
    ]
    sphingoid_found = False
    for pat in sph_base_smarts_patterns:
        sph_pattern = Chem.MolFromSmarts(pat)
        if sph_pattern and mol.HasSubstructMatch(sph_pattern):
            sphingoid_found = True
            break
    if not sphingoid_found:
        return False, "Missing sphingoid base fragment (expected carbonyl-N-two adjacent chiral carbons with OH)"
    
    # Check for a long aliphatic (fatty acid) chain.
    # We look for a chain of 8 consecutive carbons (this is a heuristic).
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long aliphatic (fatty acid) chain"
    
    # (Optional) A molecular weight check could be added if desired.
    # wt = rdMolDescriptors.CalcExactMolWt(mol)
    # if wt < 500:
    #     return False, f"Molecular weight too low ({wt:.1f} Da) for sphingomyelin"
    
    # If all checks pass, classify as sphingomyelin.
    return True, "Contains phosphorylcholine headgroup, amide-linked fatty acid, and a sphingoid base fragment"

# Example usage:
if __name__ == '__main__':
    # Example test with one sphingomyelin structure
    test_smiles = "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
    result, reason = is_sphingomyelin(test_smiles)
    print("Result:", result)
    print("Reason:", reason)