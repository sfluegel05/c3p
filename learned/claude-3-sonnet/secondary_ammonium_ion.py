"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
Definition: An organic cation obtained by protonation of any secondary amino compound; 
           major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has an NHR2+ group where R represents organic groups.

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

    # SMARTS pattern for secondary ammonium ion
    # [NH+] with exactly two single bonds to carbon-containing groups
    sec_ammonium_pattern = Chem.MolFromSmarts('[NX3H1+;!$(N=*);!$(N#*);!$([N][O,N,S,P])](-[#6])-[#6]')
    
    # Patterns to exclude
    exclude_patterns = [
        Chem.MolFromSmarts('[nH+]'), # Aromatic N+
        Chem.MolFromSmarts('[NH+]=*'), # Iminium ions
        Chem.MolFromSmarts('[NH+]#*'), # N+ triple bonds
    ]

    # Check if molecule matches the secondary ammonium pattern
    if not mol.HasSubstructMatch(sec_ammonium_pattern):
        return False, "No secondary ammonium ion found"
        
    # Check for excluding patterns
    for pattern in exclude_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains excluded nitrogen cation pattern"

    # Additional validation of each matching nitrogen
    matches = mol.GetSubstructMatches(sec_ammonium_pattern)
    for match in matches:
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Must have +1 charge
        if n_atom.GetFormalCharge() != 1:
            continue
            
        # Must not be aromatic
        if n_atom.GetIsAromatic():
            continue
            
        # Count non-hydrogen neighbors
        non_h_neighbors = [neighbor for neighbor in n_atom.GetNeighbors() 
                         if neighbor.GetAtomicNum() != 1]
        
        # Must have exactly 2 non-hydrogen neighbors
        if len(non_h_neighbors) != 2:
            continue
            
        # Verify both neighbors are carbon-containing organic groups
        if not all(neighbor.GetAtomicNum() == 6 for neighbor in non_h_neighbors):
            continue
            
        # All checks passed
        return True, "Contains secondary ammonium ion (NHR2+)"
        
    return False, "No valid secondary ammonium ion found"

__metadata__ = {
    'chemical_class': {
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['organic cation', 'ammonium ion']
    }
}