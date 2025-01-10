"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA contains a 3-hydroxy group on the fatty acid portion
    and a CoA moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a 3-hydroxy group on a fatty acyl chain
    hydroxy_acyl_pattern = Chem.MolFromSmarts("[C@@H](O)CC(=O)C")
    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "No 3-hydroxy group on a fatty acyl chain found"
    
    # Define a more comprehensive and specific SMARTS pattern for the CoA moiety
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Check additional features to minimize false positives
    # Count the number of atoms in the molecule and ensure it's within a typical range for CoAs
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 50 or num_atoms > 100:
        return False, f"Number of atoms ({num_atoms}) is not typical for a 3-hydroxy fatty acyl-CoA"
    
    return True, "Identified as a 3-hydroxy fatty acyl-CoA with CoA moiety"