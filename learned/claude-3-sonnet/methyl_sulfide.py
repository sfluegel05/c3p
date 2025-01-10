"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule contains a methyl sulfide group based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for S-CH3 pattern
    # [#16X2] = divalent sulfur
    # [CH3X4] = sp3 methyl group
    methylsulfide_pattern = Chem.MolFromSmarts("[#16X2][CH3X4]")
    matches = mol.GetSubstructMatches(methylsulfide_pattern)
    
    if not matches:
        return False, "No methylsulfide (S-CH3) group found"
    
    # Check that sulfur is not part of a higher oxidation state group
    # Look for S=O, S(=O)=O patterns that would indicate sulfoxide/sulfone
    sulfoxide_pattern = Chem.MolFromSmarts("[#16X3](=[OX1])")
    sulfone_pattern = Chem.MolFromSmarts("[#16X4](=[OX1])=[OX1]")
    
    for match in matches:
        s_atom_idx = match[0]  # Index of the sulfur atom
        s_atom = mol.GetAtomWithIdx(s_atom_idx)
        
        # Check if this sulfur is part of sulfoxide/sulfone
        if mol.HasSubstructMatch(sulfoxide_pattern) or mol.HasSubstructMatch(sulfone_pattern):
            continue
            
        # Verify sulfur has exactly two single bonds (sulfide)
        if s_atom.GetDegree() == 2 and s_atom.GetTotalValence() == 2:
            return True, "Contains methylsulfide (S-CH3) group"
            
    return False, "Sulfur is not in sulfide form (may be sulfoxide, sulfone, or other S species)"