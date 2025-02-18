"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid (contains at least one C=C double bond and a carboxylic acid group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group"
    
    # Check for at least one double bond
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds == 0:
        return False, "No double bonds present"
    
    # Optional: Check chain length (e.g., at least 8 carbons in the main chain)
    # This is a heuristic; adjust as needed
    main_chain_length = Chem.rdMolDescriptors.CalcLongestChain(mol)
    if main_chain_length < 8:
        return False, f"Main chain too short ({main_chain_length} carbons)"
    
    return True, "Contains carboxylic acid and at least one double bond"