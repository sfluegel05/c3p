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
    Determines if a molecule is a 2‑enoyl‑CoA based on its SMILES string.
    A 2‑enoyl‑CoA is defined as an unsaturated fatty acyl-CoA in which the 
    acyl (fatty acid) moiety is attached via a thioester group and contains a double bond 
    between the alpha (position 2) and beta (position 3) carbons (i.e. between the carbonyl and its adjacent carbon)
    and that the molecule contains features of CoA (e.g. an adenine nucleotide).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a 2‑enoyl‑CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for the unsaturated acyl unit:
    # We're looking for a thioester group (C(=O)S) directly followed by a carbon which is double-bonded to another carbon.
    # This allows for substitution on the alpha carbon as is common in some acyl chains.
    unsat_acyl_smarts = "C(=O)S[C]=[C]"
    unsat_acyl_pat = Chem.MolFromSmarts(unsat_acyl_smarts)
    if unsat_acyl_pat is None:
        return False, "Error building unsaturated acyl SMARTS pattern"
    
    if not mol.HasSubstructMatch(unsat_acyl_pat):
        return False, "The unsaturated acyl moiety (double bond between positions 2 and 3) was not found"
    
    # SMARTS pattern for a CoA-like moiety:
    # Rather than matching an extended pantetheine chain (which might vary due to charge/protonation states),
    # we look for the adenine ring that is part of the CoA scaffold.
    coa_smarts = "n1cnc2ncnc12"
    coa_pat = Chem.MolFromSmarts(coa_smarts)
    if coa_pat is None:
        return False, "Error building CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pat):
        return False, "CoA structural features (adenine ring) not found in the molecule"
    
    return True, "Molecule contains a thioester-linked acyl group with a double bond between positions 2 and 3 and a CoA moiety"