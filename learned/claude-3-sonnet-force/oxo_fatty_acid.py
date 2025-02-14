"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:35638 oxo fatty acid

An oxo fatty acid is any fatty acid containing at least one aldehydic or ketonic group 
in addition to the carboxylic acid group.

Typical features:
- Long aliphatic chain (6-24 carbons)
- Terminal carboxylic acid group (-C(=O)O)
- At least one aldehydic (-C(=O)H) or ketonic (-C(=O)-) group
- Molecular weight typically between 100-350 Da
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for oxo fatty acid
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("[C;$(C(=O)O)][-;!$(C=O)]" +  # Terminal carboxyl group
                                                 "[C;$(C(=O)[C,H])]" +         # Oxo group (aldehyde or ketone)
                                                 "[C;$(CC(=O)[C,H])]" +        # One carbon away from the oxo group
                                                 "[C;$(CCCC)]")                # At least 4 carbons away from the carboxyl group
    
    # Check if the molecule matches the pattern
    matches = mol.GetSubstructMatches(oxo_fatty_acid_pattern)
    
    if not matches:
        return False, "Molecule does not match the oxo fatty acid pattern"
    
    # Count carbon chain length
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check carbon chain length
    if carbon_chain_length < 6 or carbon_chain_length > 24:
        return False, f"Carbon chain length ({carbon_chain_length}) outside typical range for fatty acids (6-24)"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 350:
        return False, f"Molecular weight ({mol_wt:.2f} Da) outside typical range for oxo fatty acids (100-350 Da)"
    
    return True, "Molecule matches the structure of an oxo fatty acid"