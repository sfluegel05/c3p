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
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    
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
    
    # Look for S-CH3 pattern where S is connected to any carbon (aliphatic or aromatic)
    # [#16X2] = divalent sulfur
    # [CH3X4] = sp3 methyl group
    # [#6] = any carbon
    methylsulfide_pattern = Chem.MolFromSmarts("[#16X2]([CH3X4])[#6]")
    matches = mol.GetSubstructMatches(methylsulfide_pattern)
    
    if not matches:
        return False, "No methylsulfide (S-CH3) group found"
    
    # Check that sulfur is not part of a higher oxidation state group
    sulfoxide_pattern = Chem.MolFromSmarts("[#16X3](=[OX1])")
    sulfone_pattern = Chem.MolFromSmarts("[#16X4](=[OX1])=[OX1]")
    sulfinic_pattern = Chem.MolFromSmarts("[#16X3]([OX2H])")
    sulfonic_pattern = Chem.MolFromSmarts("[#16X4]([OX2H])(=[OX1])=[OX1]")
    
    if any(mol.HasSubstructMatch(pattern) for pattern in 
           [sulfoxide_pattern, sulfone_pattern, sulfinic_pattern, sulfonic_pattern]):
        return False, "Sulfur is in higher oxidation state"
    
    # Check if sulfur is part of a disulfide or polysulfide
    disulfide_pattern = Chem.MolFromSmarts("[SX2]-[SX2]")
    if mol.HasSubstructMatch(disulfide_pattern):
        return False, "Contains disulfide bond"
        
    # Get molecular weight
    mol_weight = sum([atom.GetMass() for atom in mol.GetAtoms()])
    
    # Filter out very small molecules (like dimethyl sulfide)
    if mol_weight < 60:
        return False, "Molecule too small to be considered a methyl sulfide compound"

    # Verify that the sulfur is connected to exactly two atoms
    for match in matches:
        s_atom_idx = match[0]  # Index of the sulfur atom
        s_atom = mol.GetAtomWithIdx(s_atom_idx)
        if s_atom.GetDegree() != 2:
            return False, "Sulfur has incorrect number of connections"
            
    # Check for thiol groups (R-SH)
    thiol_pattern = Chem.MolFromSmarts("[SX2H]")
    if mol.HasSubstructMatch(thiol_pattern):
        return False, "Contains thiol group"
        
    return True, "Contains methylsulfide (S-CH3) group with appropriate structure"