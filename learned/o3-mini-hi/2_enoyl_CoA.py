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
    A 2‑enoyl‑CoA is defined as an unsaturated fatty acyl-CoA molecule, 
    where the acyl (fatty acid) moiety contains a double bond between positions 2 and 3
    (i.e., between the α- and β-carbons counting from the carbonyl carbon) 
    and is attached to a CoA scaffold.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2‑enoyl‑CoA, False otherwise
        str: Reason for classification
    """
    # Attempt to parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS to check for the unsaturated acyl group.
    # This pattern seeks a fragment where an sp2 carbon (=CH) is directly connected 
    # to another sp2 carbon which in turn is immediately bonded to a carbonyl carbon (C(=O)S).
    # This effectively tests for a C(=O) group with the α-carbon (C2) double-bonded to the β-carbon (C3).
    unsat_acyl_smarts = "[CH]=[CH]C(=O)S"
    unsat_acyl_pattern = Chem.MolFromSmarts(unsat_acyl_smarts)
    if unsat_acyl_pattern is None:
        return False, "Error building unsaturated acyl SMARTS pattern"
    
    if not mol.HasSubstructMatch(unsat_acyl_pattern):
        return False, "The unsaturated acyl moiety (double bond between positions 2 and 3) was not found"

    # SMARTS to detect a fragment of the CoA (pantetheine) moiety.
    # Most acyl-CoA structures contain a characteristic fragment like "SCCNC(=O)CCNC(=O)"
    # that links the acyl group to the nucleotide portion.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error building CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structural features not found in the molecule"

    return True, "Molecule contains an unsaturated acyl group with a double bond between positions 2 and 3 and a CoA moiety"