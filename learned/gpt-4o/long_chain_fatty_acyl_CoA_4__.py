"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # CoA structural patterns, refined for better accuracy.
    coa_patterns = {
        'thioester': Chem.MolFromSmarts("C(=O)S"),  # Thioester linkage
        'peptide': Chem.MolFromSmarts("NC(=O)CCNC"),  # Part of peptide-like segment
        'adenine': Chem.MolFromSmarts("n1cnc2c(N)ncnc12"),  # Adenine nucleoside
        'ribose_phosphate': Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H](O)[C@@H]1OP(=O)([O-])[O-])"),  # Ribose linked with phosphate
        'dephosphorylated': Chem.MolFromSmarts("COP(=O)([O-])[O-]")  # Cashed phosphate group
    }
    
    # Check each substructure pattern in CoA
    for name, pattern in coa_patterns.items():
        if not mol.HasSubstructMatch(pattern):
            return False, f"Missing or incomplete CoA moiety: {name}"

    # Check deprotonated phosphate groups
    phos_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])[O]")
    phos_matches = mol.GetSubstructMatches(phos_pattern)
    if len(phos_matches) < 2:
        return False, "Must have at least two deprotonated phosphate groups"

    # Count carbon atoms, especially in the acyl chain
    # Consider acyl long-chain typically having at least 14 carbons
    carbon_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_count += 1

    if carbon_count < 16:
        return False, "Acyl chain is too short to be considered long-chain"

    return True, "The molecule is classified as long-chain fatty acyl-CoA(4-)"