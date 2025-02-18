"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group 
contains a double bond between positions 2 and 3.
"""

from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the 
    acyl (fatty acid) moiety is attached via a thioester group (S-C(=O)-R) and 
    the acyl chain R contains a double bond between the alpha (position 2) 
    and beta (position 3) carbons -- that is, the first two carbons of R after the carbonyl.
    In addition, the molecule should display features of CoA (such as the adenine ring).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 2-enoyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to capture the unsaturated acyl moiety attached via a thioester.
    # This pattern looks for a sulfur (S) directly bonded to a carbonyl carbon (C(=O))
    # which is bonded to a carbon that is double-bonded to another carbon.
    unsat_acyl_pat = Chem.MolFromSmarts("[S]C(=O)[C]=[C]")
    if unsat_acyl_pat is None:
        return False, "Error building unsaturated acyl SMARTS pattern"
    
    if not mol.HasSubstructMatch(unsat_acyl_pat):
        return False, "The unsaturated acyl moiety (double bond between positions 2 and 3) was not found"
    
    # SMARTS pattern for CoA: using the adenine group (n1cnc2ncnc12) as a diagnostic
    coa_pat = Chem.MolFromSmarts("n1cnc2ncnc12")
    if coa_pat is None:
        return False, "Error building CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pat):
        return False, "CoA structural features (adenine ring) not found in the molecule"
    
    return True, "Molecule contains a thioester-linked acyl group with a double bond between positions 2 and 3 and a CoA moiety"

# Example usage:
# smi = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
# result, reason = is_2_enoyl_CoA(smi)
# print(result, reason)