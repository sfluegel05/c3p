"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
#!/usr/bin/env python
"""
Classifies: Sphingomyelin d18:1 – any sphingomyelin having sphingosine as the sphingoid component.
Improved heuristic criteria:
  1. Must contain the phosphocholine head group (SMARTS: "COP(=O)([O-])OCC[N+](C)(C)C").
  2. Must contain an amide bond (SMARTS: "NC(=O)").
  3. Must NOT contain ester bonds (SMARTS: "[CX3](=O)O[CX4]").
  4. Must have a minimum overall carbon count (now relaxed to at least 35).
  5. Must have at least one isolated (nonaromatic) C=C double bond.
  6. Must contain a sphingosine-like backbone fragment, defined as a three-carbon chain where:
       - The first carbon is attached (via O) to the phosphocholine group,
       - The second bears an N-acyl substituent (amide linkage),
       - The third carries a free hydroxyl.
     We capture this using a more permissive SMARTS:
       "[C]-OP(=O)([O-])OCC[N+](C)(C)C-[C](NC(=O))-[C](O)"
     
Note: These are heuristic rules; many borderline cases exist.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines whether a molecule (given by SMILES) is a sphingomyelin d18:1.
    
    The function checks for:
      1. The phosphocholine head group.
      2. An amide bond (indicative of the N-acyl linkage).
      3. Absence of ester bonds.
      4. A sufficiently high overall carbon count (>=35).
      5. At least one isolated C=C double bond.
      6. The presence (exactly once) of a sphingosine-like backbone fragment.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets our heuristic criteria for sphingomyelin d18:1, else False.
        str: A short reason describing the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 1. Check for the phosphocholine head group.
    phos_choline_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_choline = Chem.MolFromSmarts(phos_choline_smarts)
    if not mol.HasSubstructMatch(phos_choline):
        return False, "Phosphocholine head group not found."
    
    # 2. Check for an amide bond (indicative of the N-acyl linkage).
    amide_smarts = "NC(=O)"
    amide = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide):
        return False, "Amide bond (N-acyl linkage) not found."
    
    # 3. Reject molecules that have ester bonds.
    ester_smarts = "[CX3](=O)O[CX4]"
    ester = Chem.MolFromSmarts(ester_smarts)
    if mol.HasSubstructMatch(ester):
        return False, "Ester bond detected; likely a glycerophospholipid."
    
    # 4. Check for a minimum overall carbon count (relaxed to 35).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbon atoms ({c_count}) for sphingomyelin d18:1."
    
    # 5. Check for at least one isolated (nonaromatic) C=C double bond.
    dbl_bond_smarts = "[C;!R]=[C;!R]"
    dbl_bond = Chem.MolFromSmarts(dbl_bond_smarts)
    if not mol.HasSubstructMatch(dbl_bond):
        return False, "No isolated C=C double bond found (unsaturation required)."
    
    # 6. Check for the sphingosine-like backbone fragment.
    # Using a more permissive SMARTS pattern that does not force a fixed hydrogen count or chirality.
    sphingosine_smarts = "[C]-OP(=O)([O-])OCC[N+](C)(C)C-[C](NC(=O))-[C](O)"
    sphingosine_frag = Chem.MolFromSmarts(sphingosine_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingosine_frag)
    if len(sphingo_matches) != 1:
        return False, f"Sphingosine-like backbone fragment not detected uniquely (found {len(sphingo_matches)} matches)."
    
    return True, "Molecule meets heuristic criteria for sphingomyelin d18:1 (phosphocholine head, N-acyl linkage, no ester, adequate C count, unsaturation, and sphingosine backbone)."

# Example usage for testing
if __name__ == "__main__":
    examples = [
        # Provided examples
        ("C(=C\\C/C=C\\CCCCC)\\CCCCCCCCCC(CC(=O)N[C@@H](COP(OCC[N+](C)(C)C)(=O)[O-])[C@@H](\\C=C\\CCCCCCCCCCCCC)O)O",
         "N-[(13Z,16Z)-3-hydroxydocosa-13,16-enoyl]sphingosine-1-phosphocholine"),
        ("[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCC",
         "N-pentadecanoylsphingosine-1-phosphocholine"),
        ("C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C",
         "N-docosanoylsphingosine-1-phosphocholine"),
        # Other examples can be added similarly...
    ]
    
    for smi, name in examples:
        result, reason = is_sphingomyelin_d18_1(smi)
        print("Name:", name)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 70)