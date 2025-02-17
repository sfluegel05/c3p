"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
#!/usr/bin/env python
"""
Classifies: Sphingomyelin d18:1 – any sphingomyelin having sphingosine as the sphingoid component.
Heuristic criteria:
  1. Must contain the phosphocholine head group (SMARTS: "COP(=O)([O-])OCC[N+](C)(C)C").
  2. Must contain an amide bond (SMARTS: "NC(=O)").
  3. Must not have any ester bonds (SMARTS: "[CX3](=O)O[CX4]") that would indicate diacyl glycerophospholipids.
  4. Must possess a minimum overall number of carbons (here at least 35).
  5. Must possess at least one isolated C=C double bond (nonaromatic), as a proxy for unsaturation in the sphingosine backbone.
Note: These criteria are heuristic and do not capture every nuance.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 is defined as any sphingomyelin having sphingosine as the sphingoid component.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as sphingomyelin d18:1, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # 1. Check for phosphocholine head group.
    # The SMARTS looks for: –O–P(=O)([O-])–O–CC[N+](C)(C)C.
    phos_choline_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_choline = Chem.MolFromSmarts(phos_choline_smarts)
    if not mol.HasSubstructMatch(phos_choline):
        return False, "Phosphocholine head group not found."

    # 2. Check for an amide bond.
    # The amide bond should appear as N followed directly by a carbonyl group.
    amide_smarts = "NC(=O)"
    amide = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide):
        return False, "Amide bond (fatty acyl linkage) not found."

    # 3. Exclude molecules containing ester bonds.
    # Diacyl glycerophosphocholines have ester bonds (i.e. C(=O)O–C) connecting acyl chains to glycerol.
    # We use the SMARTS "[CX3](=O)O[CX4]". (Note: The ester linkage in the head group is a phosphate ester and does not match this pattern.)
    ester_smarts = "[CX3](=O)O[CX4]"
    ester = Chem.MolFromSmarts(ester_smarts)
    if mol.HasSubstructMatch(ester):
        return False, "Ester bond detected; molecule likely is a diacyl glycerophosphocholine rather than sphingomyelin d18:1."

    # 4. Check for minimum overall carbon count.
    # Sphingomyelins consist of an 18-carbon sphingosine base, a fatty acyl chain, and the head group.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbon atoms ({c_count}) for a sphingomyelin d18:1."

    # 5. Check for at least one isolated (nonaromatic) C=C double bond.
    # This is a proxy for unsaturation in the sphingosine backbone.
    dbl_bond_smarts = "[C;!R]=[C;!R]"
    dbl_bond = Chem.MolFromSmarts(dbl_bond_smarts)
    if not mol.HasSubstructMatch(dbl_bond):
        return False, "No isolated C=C double bond found."

    # If all criteria are met, we're classifying the molecule as sphingomyelin d18:1.
    return True, "Molecule contains phosphocholine head group, an amide linkage, no extraneous ester bonds, and satisfies other heuristic criteria for sphingomyelin d18:1."

# Example usage:
if __name__ == "__main__":
    # List of example SMILES for molecules that should be classified as sphingomyelin d18:1.
    examples = [
        "C(=C\\C/C=C\\CCCCC)\\CCCCCCCCCC(CC(=O)N[C@@H](COP(OCC[N+](C)(C)C)(=O)[O-])[C@@H](\\C=C\\CCCCCCCCCCCCC)O)O",
        "[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCC",
        "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    ]
    
    for smi in examples:
        result, reason = is_sphingomyelin_d18_1(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 50)