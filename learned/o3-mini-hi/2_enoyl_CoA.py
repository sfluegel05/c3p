"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
This program checks whether the molecule has:
  1. A thioester linkage (S-C(=O)-...) where the carbonyl carbon is attached to an acyl chain that has a double bond between the first and second carbons (positions 2 and 3 relative to the carbonyl).
  2. A CoA moiety, at least as indicated by the presence of an adenine ring.
The code removes stereochemistry to make the substructure search more robust.
"""

from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the 
    acyl group is attached via a thioester linkage (S-C(=O)-R), and the acyl chain
    R bears a double bond between the alpha and beta carbons (positions 2 and 3 relative to the carbonyl).
    In addition, the molecule should contain CoA structural features (e.g. an adenine ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 2-enoyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemistry to avoid mismatch due to chiral annotations.
    Chem.RemoveStereochemistry(mol)
    
    # Define a SMARTS pattern for the unsaturated acyl (thioester) group.
    # This pattern looks for a sulfur atom bonded to a carbonyl carbon,
    # which in turn is bonded to a carbon that is double-bonded to another carbon.
    unsat_acyl_pat = Chem.MolFromSmarts("[S]C(=O)[C]=[C]")
    if unsat_acyl_pat is None:
        return False, "Error building unsaturated acyl SMARTS pattern"

    if not mol.HasSubstructMatch(unsat_acyl_pat):
        return False, "The unsaturated acyl moiety with a double bond between positions 2 and 3 was not found"
    
    # Define a pattern for the CoA moiety.
    # We use the adenine fragment as a diagnostic. The canonical SMILES for adenine is:
    # "c1nc2c(n1)nc[nH]2". Using this as a query tends to be more flexible.
    adenine = Chem.MolFromSmiles("c1nc2c(n1)nc[nH]2")
    if adenine is None:
        return False, "Error building adenine substructure pattern for CoA"
    
    if not mol.HasSubstructMatch(adenine):
        return False, "CoA structural features (adenine ring) not found in the molecule"
    
    return True, "Molecule contains a thioester-linked acyl group with a double bond between positions 2 and 3 and displays CoA features (adenine ring)"

# Example usage:
if __name__ == "__main__":
    # Example SMILES for (E)-hexadec-2-enoyl-CoA; you can add more tests.
    smi = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_2_enoyl_CoA(smi)
    print(result, reason)