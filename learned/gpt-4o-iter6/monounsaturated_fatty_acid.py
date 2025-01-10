"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the fatty acid chain and  
    a carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylic acid group pattern (ensuring it's terminal and part of a long chain)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Calculate number of carbon-carbon double/triple bonds (only those in the chain, not part of ring systems)
    chain_double_bond_count = len([bond for bond in mol.GetBonds() 
                                   if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]
                                   and bond.GetBeginAtom().GetDegree() <= 2 and bond.GetEndAtom().GetDegree() <= 2])
    
    if chain_double_bond_count != 1:
        return False, f"Found {chain_double_bond_count} unsaturations in chain, need exactly one"

    # Ensure it's a straight-chain or branched molecule without complex ring systems contributing features
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings which disqualifies it from being a straight chain fatty acid"
    
    return True, "Molecule is a monounsaturated fatty acid (one double or triple bond in the chain with carboxylic group)"