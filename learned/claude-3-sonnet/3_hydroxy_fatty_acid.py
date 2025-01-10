"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for specific 3-hydroxy fatty acid pattern
    # Pattern matches: OH-CH(R)-CH2-CH2-COOH or OH-CH2-CH2-CH2-COOH
    beta_hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]-[CX4][CX4][CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No hydroxy group at 3-position"

    # Exclude peptides by checking for absence of amino groups near the acid
    peptide_pattern = Chem.MolFromSmarts("[NX3,NX4]-[CX4]-[CX3](=[OX1])[OX2H1]")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Appears to be a peptide"
        
    # Exclude glucuronides and other sugar derivatives
    sugar_pattern = Chem.MolFromSmarts("[OX2]-[CX4]-[OX2]-[CX4]-[OX2]-[CX4]")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moiety"

    # Count carbons to ensure it's a fatty acid (minimum 5 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Count oxygens (should typically have 3-4 oxygens: COOH + OH)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms for 3-hydroxy fatty acid"
    if oxygen_count > 8:  # Allow for some additional OH groups but not too many
        return False, "Too many oxygen atoms for typical 3-hydroxy fatty acid"

    # Verify the molecule has a predominant carbon chain
    # Count carbons with more than 2 carbon neighbors
    highly_branched_carbons = sum(1 for atom in mol.GetAtoms() 
                                if atom.GetAtomicNum() == 6 and 
                                len([n for n in atom.GetNeighbors() 
                                    if n.GetAtomicNum() == 6]) > 2)
    if highly_branched_carbons > carbon_count / 5:  # Allow some branching but not too much
        return False, "Carbon chain too branched for typical fatty acid"

    # Check for aromatic rings which shouldn't be present in fatty acids
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Count rotatable bonds to verify chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too rigid to be a fatty acid"

    return True, "Contains 3-hydroxy group and appropriate fatty acid chain characteristics"