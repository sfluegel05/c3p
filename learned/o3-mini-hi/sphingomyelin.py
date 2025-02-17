"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: Sphingomyelin – any phospholipid where the amino group of a sphingoid base 
is in amide linkage with a fatty acid, and the terminal hydroxy group is esterified 
to phosphorylcholine.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelins must have an amide bond linking a fatty acid to a sphingoid base
    and the terminal hydroxyl of the base is esterified to a phosphorylcholine headgroup.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as sphingomyelin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the phosphorylcholine headgroup.
    # This pattern looks for a –O–P(=O)([O-])–OCC[N+](C)(C)C fragment.
    phos_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phos_pattern):
        return False, "Missing phosphorylcholine headgroup"
    
    # Define a SMARTS pattern for a basic amide linkage (fatty acid part linked via an amide).
    # We look for an N atom bonded to a carbonyl group.
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide bond linking a fatty acid to the sphingoid base"
    
    # Optional: Check that the molecule has a sufficiently high molecular weight 
    # (sphingomyelins are typically >600 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for sphingomyelin"
    
    # Optional: Count the total number of carbon atoms (should have many aliphatic carbons).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:
        return False, "Insufficient number of carbons for a long-chain sphingomyelin structure"
        
    # If all checks pass, we classify the molecule as sphingomyelin.
    return True, "Contains amide-linked fatty acid and phosphorylcholine ester as required for sphingomyelin"

# Example usage:
if __name__ == '__main__':
    # Test with one example SMILES (N-hexadecanoyl-15-methylhexadecasphing-4-enine-1-phosphocholine)
    test_smiles = "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCC(C)C"
    result, reason = is_sphingomyelin(test_smiles)
    print("Result:", result)
    print("Reason:", reason)