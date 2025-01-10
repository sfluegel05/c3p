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
    
    # Look for S-CH3 pattern where S is connected to another carbon
    # [#16X2] = divalent sulfur
    # [CH3X4] = sp3 methyl group
    # [#6X4] = sp3 carbon (aliphatic)
    methylsulfide_pattern = Chem.MolFromSmarts("[#16X2]([CH3X4])[#6X4]")
    matches = mol.GetSubstructMatches(methylsulfide_pattern)
    
    if not matches:
        return False, "No aliphatic methylsulfide (S-CH3) group found"
    
    # Check that sulfur is not part of a higher oxidation state group
    sulfoxide_pattern = Chem.MolFromSmarts("[#16X3](=[OX1])")
    sulfone_pattern = Chem.MolFromSmarts("[#16X4](=[OX1])=[OX1]")
    
    if mol.HasSubstructMatch(sulfoxide_pattern) or mol.HasSubstructMatch(sulfone_pattern):
        return False, "Sulfur is in higher oxidation state (sulfoxide/sulfone)"
    
    # Get molecular weight and atom count
    mol_weight = sum([atom.GetMass() for atom in mol.GetAtoms()])
    atom_count = mol.GetNumAtoms()
    
    # Filter out very small molecules (like dimethyl sulfide)
    if mol_weight < 80 or atom_count < 5:
        return False, "Molecule too small to be considered a methyl sulfide compound"
        
    # Filter out peptides containing methionine
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
        if peptide_bonds >= 1:
            return False, "Peptide containing methionine, not a methyl sulfide compound"
    
    # Check if sulfur is part of a ring
    ri = mol.GetRingInfo()
    for match in matches:
        s_atom_idx = match[0]  # Index of the sulfur atom
        if ri.IsAtomInRingOfSize(s_atom_idx, 3) or ri.IsAtomInRingOfSize(s_atom_idx, 4):
            return False, "Sulfur is part of a small ring system"
            
    return True, "Contains aliphatic methylsulfide (S-CH3) group with appropriate molecular complexity"