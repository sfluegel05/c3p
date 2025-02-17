"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
This program checks whether the molecule has:
  1. A thioester linkage (S-C(=O)-...) where the carbonyl carbon is connected to an acyl chain that has a double bond between the alpha and beta carbons.
  2. A CoA moiety, at least as indicated by the presence of an adenine ring.
The code removes stereochemistry to make the substructure search more robust.
"""

from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the 
    acyl group is attached via a thioester linkage (S-C(=O)-R), and the acyl chain
    R bears a double bond between the alpha and beta carbons (positions 2 and 3 relative to the carbonyl).
    Additionally, the molecule should have features of CoA (i.e. an adenine ring).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a 2-enoyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemistry to avoid mismatches due to chiral annotations.
    Chem.RemoveStereochemistry(mol)
    
    # Define a SMARTS pattern for the unsaturated acyl group.
    # We look for a sulfur atom attached to a carbonyl carbon which in turn connects to a carbon-carbon double bond.
    unsat_acyl_pat = Chem.MolFromSmarts("[S]C(=O)[C]=[C]")
    if unsat_acyl_pat is None:
        return False, "Error building unsaturated acyl SMARTS pattern"
    
    if not mol.HasSubstructMatch(unsat_acyl_pat):
        return False, ("The unsaturated acyl moiety with a thioester linkage and "
                       "a double bond between positions 2 and 3 was not found")
    
    # Define a SMARTS pattern for the adenine ring, as a proxy for CoA's adenine moiety.
    # Using a variant that has been shown to build reliably:
    adenine_pat = Chem.MolFromSmarts("n1cnc2ncnc12")
    if adenine_pat is None:
        return False, "Error building adenine substructure pattern for CoA"
    
    if not mol.HasSubstructMatch(adenine_pat):
        return False, "CoA features (adenine ring) not found in the molecule"
    
    return True, ("Molecule contains a thioester-linked acyl group with a double bond between "
                  "positions 2 and 3 and displays CoA features (adenine ring)")

# Example usage:
if __name__ == "__main__":
    # Example SMILES for (E)-hexadec-2-enoyl-CoA.
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_2_enoyl_CoA(test_smiles)
    print(result, reason)