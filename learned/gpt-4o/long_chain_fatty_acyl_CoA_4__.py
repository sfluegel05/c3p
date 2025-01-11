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
    
    # Identify CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("C[C@H](O)[C@](C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n2cnc3c(ncnc23)N)N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA moiety"
    
    # Identify deprotonated phosphate groups
    phos_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])[O]")
    phos_matches = mol.GetSubstructMatches(phos_pattern)
    if len(phos_matches) < 2:
        return False, "Must have deprotonated phosphate groups"

    # Check for long fatty acyl chain (aliphatic connected to the thioester)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage for acyl connection"

    # Estimate the number of carbons in the acyl chain
    acyl_chain = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC(=O)")
    if not mol.HasSubstructMatch(acyl_chain):
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count < 16:
            return False, "Acyl chain is too short to be long-chain"

    return True, "Contains long-chain fatty acyl-CoA(4-) components"