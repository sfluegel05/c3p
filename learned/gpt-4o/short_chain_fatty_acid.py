"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with fewer than 6 carbon atoms,
    and without any non-hydrocarbon substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carboxylic acid group: -C(=O)O
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count >= 6:
        return False, f"Contains {c_count} carbons, expected fewer than 6"
    
    # Check for non-hydrocarbon substituents (except allowed ones like OH)
    allowed_substituents = {8}  # 8 is oxygen, like in -OH groups. Ensure no other heteroatoms.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in {1, 6} and atomic_num not in allowed_substituents:
            return False, f"Contains non-hydrocarbon substituent: {atom.GetSymbol()}"
            
    # Ensure it's aliphatic: no aromatic rings
    if mol.GetAromaticityModel() != rdchem.AromaticityModel.NonAromatic:
        return False, "Contains aromatic structures, expected aliphatic"
    
    return True, "Valid short-chain fatty acid with carboxylic acid group and <6 carbon atoms"

# Example call 
# This should return (True, "Valid short-chain fatty acid with carboxylic acid group and <6 carbon atoms")
# print(is_short_chain_fatty_acid("CCCC(O)=O"))  # Represents butyric acid