"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: Sphingomyelin d18:1 – any sphingomyelin having sphingosine as the sphingoid component.
Heuristic criteria used here:
  1. The molecule must contain the phosphocholine headgroup (SMARTS: "COP(=O)([O-])OCC[N+](C)(C)C").
  2. It must contain an amide bond (SMARTS: "NC(=O)") which connects a fatty acyl chain to the sphingoid base.
  3. The overall number of carbons is expected to be high (here we require at least 35 carbons).
  4. There is at least one isolated C=C double bond present.
Note: These criteria are heuristic and may not cover all cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) qualifies as a sphingomyelin d18:1.
    (That is, it must contain a phosphocholine head group and an amide linkage and by assumption
    have a sphingosine d18:1 backbone.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as sphingomyelin d18:1, False otherwise.
        str: Reason for the classification.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Check for phosphocholine headgroup.
    # This pattern looks for the characteristic fragment: -O-C-P(=O)([O-])-O-CC[N+](C)(C)C
    phos_choline_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_choline = Chem.MolFromSmarts(phos_choline_smarts)
    if not mol.HasSubstructMatch(phos_choline):
        return False, "Phosphocholine head group not found."
    
    # Check for an amide bond (typical for sphingomyelin, connecting a fatty acyl chain to the sphingoid base)
    amide_smarts = "N C(=O)"
    amide = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide):
        return False, "Amide bond (fatty acyl linkage) not found."
    
    # Check that the molecule has a minimum number of carbon atoms.
    # Sphingomyelins consist of a sphingoid base (~18 carbons), a fatty acyl chain (often >14 carbons),
    # and the head group.  We require at least 35 carbons in total.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbon atoms ({c_count}) for a sphingomyelin d18:1."
    
    # Check for at least one double bond (C=C). (This is a weak proxy for unsaturation in the sphingosine backbone.)
    double_bond_smarts = "[C;!R]=[C;!R]"  # non-ring C=C (helps avoid aromatic bonds)
    dbl_bond_pattern = Chem.MolFromSmarts(double_bond_smarts)
    if not mol.HasSubstructMatch(dbl_bond_pattern):
        return False, "No isolated C=C double bond found."
        
    # We could further try to extract the specific sphingoid chain and count its carbons and double bonds,
    # but that requires additional graph analysis. Here we assume that the presence of the phosphocholine head group,
    # a fatty acyl amide linkage, sufficient overall carbon count, and at least one double bond are consistent
    # with a sphingomyelin d18:1.
    
    return True, "Molecule contains phosphocholine head group, an amide linkage, and satisfies heuristic criteria for sphingomyelin d18:1."

# Example usage:
if __name__ == "__main__":
    # List of example SMILES that (by our criteria) should be accepted as sphingomyelin d18:1.
    examples = [
        "C(=C\\C/C=C\\CCCCC)\\CCCCCCCCCC(CC(=O)N[C@@H](COP(OCC[N+](C)(C)C)(=O)[O-])[C@@H](\\C=C\\CCCCCCCCCCCCC)O)O",
        "[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCC",
        "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
        # (Add additional examples if needed)
    ]
    
    for smi in examples:
        result, reason = is_sphingomyelin_d18_1(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 40)