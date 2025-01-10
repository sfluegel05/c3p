"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule fits the class, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the CoA moiety - more flexible pattern
    coa_pattern = Chem.MolFromSmarts("NC1=NC=CN=C1N2C=C(C(=N2)N)[C@H]3O[C@H]([C@H](O)[C@@H]3OP([O-])([O-])=O)COP([O-])([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing or misidentified CoA moiety"
    
    # Identify deprotonated phosphate groups
    phos_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])[O]")
    phos_matches = mol.GetSubstructMatches(phos_pattern)
    if len(phos_matches) < 2:
        return False, "Must have deprotonated phosphate groups"

    # Check for thioester linkage as part of acyl-CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage for fatty acyl"

    # Estimate the number of carbons in the acyl chain (flexible for variability in chain length)
    acyl_chain = Chem.MolFromSmarts("C(=O)SC(C)C")
    if not mol.HasSubstructMatch(acyl_chain):
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count < 16:
            return False, "Acyl chain is too short to be long-chain"
    
    return True, "Contains long-chain fatty acyl-CoA(4-) components"