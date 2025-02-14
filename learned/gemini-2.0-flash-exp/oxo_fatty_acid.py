"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is a fatty acid (long chain carboxylic acid) containing at least one
    additional ketone or aldehyde group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
       return False, "No carboxylic acid group found"
    
    # Check for ketone or aldehyde group, not part of the carboxylic acid or a ring
    # Ketone: A carbonyl carbon connected to two other carbon atoms, that are also not carbonyls or part of the carboxylic acid.
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])([CX4;!$([CX3]=O)]X4)") 
    # Aldehyde: a carbon with one hydrogen and double bond to oxygen (not in a ring)
    aldehyde_pattern = Chem.MolFromSmarts("[CH1X3](=O)!@*")
    
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    if not (has_ketone or has_aldehyde):
       return False, "No ketone or aldehyde group found"

    #Check for linear carbon chain (at least 4 carbons) 
    # Fatty acids are typically linear, but this is not always true so we search for chains of at least 4 carbons directly connected (but allow branches).
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4;!$([CX3]=O)]~[CX4;!$([CX3]=O)]~[CX4;!$([CX3]=O)]~[CX4;!$([CX3]=O)]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
       return False, "No long carbon chain found"

   # Check if the molecule has any rings (fatty acids are generally not cyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a fatty acid"
       
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for a long carbon chain (at least 4 carbons)    
    if c_count < 4:
         return False, "Too few carbons for a fatty acid"

    # Count rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2: # Reduced rotatable bond threshold
        return False, "Too few rotatable bonds for a fatty acid"

    # Molecular weight check for fatty acid (typically > 100)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120: # Reduced minimal MW
      return False, "Molecular weight too low for fatty acid"

    return True, "Contains a carboxylic acid, a ketone or aldehyde, and a long carbon chain"