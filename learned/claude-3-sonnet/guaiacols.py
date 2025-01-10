"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule contains a guaiacol moiety based on its SMILES string.
    A guaiacol is a phenol with a methoxy group in the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains guaiacol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for guaiacol core:
    # - aromatic carbon ring
    # - hydroxyl group (-OH)
    # - methoxy group (-OCH3) in ortho position
    guaiacol_pattern = Chem.MolFromSmarts('c1c(O)c(OC)cccc1')
    alt_pattern = Chem.MolFromSmarts('c1c(OC)c(O)cccc1')  # Alternative arrangement
    
    # Check for matches
    if mol.HasSubstructMatch(guaiacol_pattern) or mol.HasSubstructMatch(alt_pattern):
        # Count matches to ensure aromatic ring
        aromatic_carbons = sum(1 for atom in mol.GetAtoms() 
                             if atom.GetIsAromatic() and atom.GetAtomicNum() == 6)
        if aromatic_carbons < 6:
            return False, "Contains required groups but not on aromatic ring"
            
        # Verify methoxy group is truly -OCH3
        methoxy_pattern = Chem.MolFromSmarts('[OH0X2]-[CH3X4]')
        if not mol.HasSubstructMatch(methoxy_pattern):
            return False, "Contains ether group but not methoxy (-OCH3)"
            
        # Verify hydroxyl group
        hydroxyl_pattern = Chem.MolFromSmarts('[OH1X2]')
        if not mol.HasSubstructMatch(hydroxyl_pattern):
            return False, "Contains oxygen but not hydroxyl (-OH)"
            
        return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"
        
    return False, "Does not contain guaiacol pattern"