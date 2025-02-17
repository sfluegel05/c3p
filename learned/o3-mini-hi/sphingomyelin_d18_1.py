"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
#!/usr/bin/env python
"""
Classifies: Sphingomyelin d18:1 – any sphingomyelin having sphingosine as the sphingoid component.

Heuristic criteria applied in this version:
  1. Must contain the phosphocholine head group (SMARTS: "COP(=O)([O-])OCC[N+](C)(C)C").
  2. Must contain an amide bond (SMARTS: "NC(=O)").
  3. Must NOT contain ester bonds indicative of glycerophospholipids (SMARTS: "[CX3](=O)O[CX4]").
  4. Must have a minimum overall number of carbon atoms (here, at least 35).
  5. Must have at least one isolated (i.e. nonaromatic) C=C double bond.
  6. Must contain a sphingosine-like backbone – a fragment where a chiral carbon bears the phosphocholine head group,
     is directly attached to a second chiral carbon that carries the amide substituent (i.e. the fatty acyl chain),
     followed by a third chiral carbon bearing an –OH group.
     
Note: These are heuristic rules and may not capture every nuance of sphingomyelin d18:1.
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
    # Try to parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 1. Check for the phosphocholine head group.
    # SMARTS that matches –O–P(=O)([O-])–O–CC[N+](C)(C)C
    phos_choline_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_choline = Chem.MolFromSmarts(phos_choline_smarts)
    if not mol.HasSubstructMatch(phos_choline):
        return False, "Phosphocholine head group not found."
    
    # 2. Check for an amide bond (indicative of the N-acyl linkage).
    amide_smarts = "NC(=O)"
    amide = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide):
        return False, "Amide bond (N-acyl linkage) not found."
    
    # 3. Reject molecules that have ester bonds (i.e. acyl chains attached to glycerol).
    ester_smarts = "[CX3](=O)O[CX4]"
    ester = Chem.MolFromSmarts(ester_smarts)
    if mol.HasSubstructMatch(ester):
        return False, "Ester bond detected; likely a glycerophospholipid rather than sphingomyelin d18:1."
    
    # 4. Check for a minimum overall carbon count.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbon atoms ({c_count}) for a typical sphingomyelin d18:1."
    
    # 5. Check for at least one isolated (nonaromatic) C=C double bond.
    dbl_bond_smarts = "[C;!R]=[C;!R]"
    dbl_bond = Chem.MolFromSmarts(dbl_bond_smarts)
    if not mol.HasSubstructMatch(dbl_bond):
        return False, "No isolated C=C double bond found (required for unsaturation in sphingosine backbone)."
    
    # 6. Check for the sphingosine backbone.
    # Heuristic SMARTS: A chiral carbon bearing the phosphocholine head group
    # linked to a chiral carbon that bears an amide substituent, followed by a chiral carbon with an –OH.
    sphingosine_smarts = "[C@@H](OP(=O)([O-])OCC[N+](C)(C)C)-[C@@H](NC(=O))-[C@@H](O)"
    sphingosine_frag = Chem.MolFromSmarts(sphingosine_smarts)
    if not mol.HasSubstructMatch(sphingosine_frag):
        return False, "Sphingosine backbone (with phosphocholine attachment, N-acyl linkage and hydroxyl) not detected."
    
    # If all heuristic criteria are met, classify as sphingomyelin d18:1.
    return True, "Molecule contains the phosphocholine head group, N-acyl linkage, no extraneous ester bonds, adequate carbon count, unsaturation, and a sphingosine-like backbone."

# Example usage:
if __name__ == "__main__":
    # Test examples that should be classified as sphingomyelin d18:1.
    examples = [
        # N-[(13Z,16Z)-3-hydroxydocosa-13,16-enoyl]sphingosine-1-phosphocholine
        "C(=C\\C/C=C\\CCCCC)\\CCCCCCCCCC(CC(=O)N[C@@H](COP(OCC[N+](C)(C)C)(=O)[O-])[C@@H](\\C=C\\CCCCCCCCCCCCC)O)O",
        # N-pentadecanoylsphingosine-1-phosphocholine
        "[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCC",
        # N-docosanoylsphingosine-1-phosphocholine
        "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C"
    ]
    
    for smi in examples:
        result, reason = is_sphingomyelin_d18_1(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)