"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin – any phospholipid where the amino group of a sphingoid base 
is in amide linkage with a fatty acid, and the terminal hydroxyl group is esterified 
to a phosphorylcholine headgroup.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelins must have an amide bond linking a fatty acid to a sphingoid base,
    the sphingoid base must have an amino group linked to a hydroxylated chiral center,
    and the terminal hydroxyl group should be esterified to a phosphorylcholine headgroup.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as sphingomyelin, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphorylcholine headgroup.
    # We allow two variants: one with explicit negative charge on an oxygen,
    # and one without the charge specification.
    phos_smarts_list = [
        "P(=O)([O-])OCC[N+](C)(C)C",  # with explicit negative charge on one oxygen
        "P(=O)(O)OCC[N+](C)(C)C"       # without explicit charge information
    ]
    phos_found = False
    for pattern in phos_smarts_list:
        phos_pattern = Chem.MolFromSmarts(pattern)
        if phos_pattern and mol.HasSubstructMatch(phos_pattern):
            phos_found = True
            break
    if not phos_found:
        return False, "Missing phosphorylcholine headgroup"

    # Check for an amide bond linking a fatty acid to the sphingoid base.
    # This pattern finds an N atom linked to a carbonyl.
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage for fatty acid attachment"

    # Check for a sphingoid base: the amino group must be attached to a chiral carbon bearing an OH.
    # We search for either of the two possible chiral notations.
    sph_base_pattern1 = Chem.MolFromSmarts("N[C@H](O)")
    sph_base_pattern2 = Chem.MolFromSmarts("N[C@@H](O)")
    if not (mol.HasSubstructMatch(sph_base_pattern1) or mol.HasSubstructMatch(sph_base_pattern2)):
        return False, "Missing sphingoid base (no N attached to a chiral center with OH)"

    # Check for a long aliphatic (fatty acyl) chain.
    # We require at least 8 consecutive carbons to indicate a long fatty acid chain.
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long aliphatic (fatty acid) chain"
    
    # (Optional) The molecular weight threshold is not enforced to allow natural variance within sphingomyelin.
    # mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # if mol_wt < 500:
    #     return False, f"Molecular weight too low ({mol_wt:.1f} Da) for sphingomyelin"
    
    # If all checks pass, classify as sphingomyelin.
    return True, "Contains amide-linked fatty acid, a sphingoid base with hydroxyl, and phosphorylcholine headgroup"

# Example usage:
if __name__ == '__main__':
    # Example test with one sphingomyelin structure
    test_smiles = "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
    result, reason = is_sphingomyelin(test_smiles)
    print("Result:", result)
    print("Reason:", reason)