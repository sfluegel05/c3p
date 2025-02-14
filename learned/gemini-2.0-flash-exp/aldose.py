"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 3:
        return False, "Aldoses must contain at least 3 carbon atoms"


    # Check for aldehyde group using SMARTS
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        
        # Verify that the carbon chain has 2 carbons, each with an -OH
        carbon_chain_pattern = Chem.MolFromSmarts("[CX4H][OX2H]~[CX4H][OX2H]")
        if mol.HasSubstructMatch(carbon_chain_pattern):
             return True, "Has an aldehyde group and at least two carbons with OH groups"
        
    else:
    # Check for intramolecular hemiacetal: check for rings containing -O-CH(-O)-
         hemiacetal_pattern = Chem.MolFromSmarts("[OX2]-[CX4H](-[OX2H])")
         if mol.HasSubstructMatch(hemiacetal_pattern):
            
            #Verify the presence of at least 2 OH groups attached to carbons
            polyhydroxy_pattern = Chem.MolFromSmarts("[CX4H][OX2H]~[CX4H][OX2H]")
            if mol.HasSubstructMatch(polyhydroxy_pattern):
                return True, "Intramolecular hemiacetal structure detected, and contains two carbons with OH groups"
    
    return False, "Does not match aldose criteria"