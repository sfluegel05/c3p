"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Simplified Coenzyme A moiety pattern
    coa_pattern = Chem.MolFromSmarts("NCC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A component found"

    # CoA also has a ribose-phosphate structure, ensuring this might help.
    coa_ribose_pattern = Chem.MolFromSmarts("OCC1OC(C(O)C1O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_ribose_pattern):
        return False, "No ribose-phosphate structure indicating CoA found"
    
    # Check for medium-chain length (6 to 12 carbons) in complete hydrocarbon chain
    carbon_chain_length = get_fatty_acid_chain_length(mol)
    if carbon_chain_length < 6 or carbon_chain_length > 12:
        return False, f"Aliphatic chain length of {carbon_chain_length} not within medium-chain range (6-12 carbons)"

    return True, "Molecule is a medium-chain fatty acyl-CoA with CoA moiety and appropriate chain length"

def get_fatty_acid_chain_length(mol):
    """
    Utility to count the longest continuous carbon chain that could represent fatty acid part.
    """
    chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 2:  # Typically part of a carbon chain
            path_lengths = Chem.rdmolops.GetShortestPaths(mol, atom.GetIdx())
            chains.extend(max(length for length in path_lengths))
    return max(chains) if chains else 0

# Example test
example_smiles = "CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
print(is_medium_chain_fatty_acyl_CoA(example_smiles))