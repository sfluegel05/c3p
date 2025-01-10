"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a volatile organic compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count types of atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])  # F, Cl, Br, I
    
    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)
    
    # Determine if it's a VOC
    # Basic check: must have more carbon than any other element and not exceed a certain molecular weight
    # Set an arbitrary upper limit to 300 g/mol to assume volatility; larger molecules tend to be less volatile
    if carbon_count >= 1 and carbon_count > oxygen_count + nitrogen_count + halogen_count and mol_wt <= 300:
        return True, "Likely a volatile organic compound based on moderate molecular weight and composition"
    else:
        return False, "Does not match basic characteristics of a volatile organic compound"