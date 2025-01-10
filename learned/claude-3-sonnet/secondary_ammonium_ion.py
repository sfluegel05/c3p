"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has an NH2+ group with exactly two carbons attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms
    potential_ammonium = False
    reason = "No secondary ammonium ion found"
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            # Check formal charge
            if atom.GetFormalCharge() != 1:
                continue
                
            # Count attached hydrogens
            n_hydrogens = atom.GetTotalNumHs()
            
            # Count attached carbons
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() 
                                if neighbor.GetAtomicNum() == 6)
            
            # Secondary ammonium should have:
            # - Exactly 2 hydrogens (NH2+)
            # - Exactly 2 carbon neighbors
            # - Formal charge of +1
            if n_hydrogens == 2 and carbon_neighbors == 2:
                return True, "Contains NH2+ group with exactly two carbons attached"
            
            # Help identify why classification failed
            if n_hydrogens > 2:
                reason = "Primary ammonium ion (NH3+)"
            elif n_hydrogens == 1 and carbon_neighbors == 3:
                reason = "Tertiary ammonium ion (NHR3+)"
            elif n_hydrogens == 0 and carbon_neighbors == 4:
                reason = "Quaternary ammonium ion (NR4+)"
            elif carbon_neighbors < 2:
                reason = f"Only {carbon_neighbors} carbon(s) attached to NH2+ group"
            elif n_hydrogens < 2:
                reason = f"Only {n_hydrogens} hydrogen(s) on charged nitrogen"
                
    return False, reason

__metadata__ = {
    'chemical_class': {
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['organic cation', 'ammonium ion']
    }
}