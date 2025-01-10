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

    # Expanded CoA moiety pattern. CoA is a complex structure with multiple amide, phosphate, and ribose components.
    coa_patterns = [
        Chem.MolFromSmarts("SCCNC"),  # Part involving thioester linkage
        Chem.MolFromSmarts("NC(=O)CCNC"),  # Peptide-like moiety of CoA
        Chem.MolFromSmarts("C(C)(C)COP(=O)([O-])[O-]"),  # Dephosphorylated phosphate group and ribose section
        Chem.MolFromSmarts("O[C@@H]1[C@@H](O)C[C@H]1OP([O-])([O-])=O"),  # Ribose connected to phosphate
        Chem.MolFromSmarts("n1cnc2c(N)ncnc12")  # Adenine nucleoside
    ]
    
    # Check each pattern in CoA
    for pattern in coa_patterns:
        if not mol.HasSubstructMatch(pattern):
            return False, "Missing or incomplete CoA moiety"

    # Identify deprotonated phosphate groups
    phos_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])[O]")
    phos_matches = mol.GetSubstructMatches(phos_pattern)
    if len(phos_matches) < 2:
        return False, "Must have deprotonated phosphate groups"

    # Check for thioester linkage as part of acyl-CoA; it might have variations.
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage for fatty acyl"
    
    # Estimate the number of carbons in the acyl chain - flexible for variability in chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
        return False, "Acyl chain is too short to be long-chain"

    return True, "Contains long-chain fatty acyl-CoA(4-) components"