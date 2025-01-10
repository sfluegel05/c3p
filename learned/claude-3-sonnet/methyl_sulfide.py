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
    
    # Look for S-CH3 pattern
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

    # Count peptide bonds
    peptide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])[#6]")
    peptide_matches = len(mol.GetSubstructMatches(peptide_pattern))
    
    # If molecule has multiple peptide bonds, it's likely a peptide
    if peptide_matches >= 2:
        return False, "Appears to be a peptide"
    
    # Check for peptide bonds near the sulfur
    s_ch3_near_peptide = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])[#6][#6][#6][SX2][CH3]")
    if mol.HasSubstructMatch(s_ch3_near_peptide):
        return False, "Methylsulfide group is part of a peptide structure"
    
    # Allow specific valid patterns like thiohydroximates
    thiohydroximate_pattern = Chem.MolFromSmarts("[#6][SX2][CH3].[#6]=[N][OH]")
    if mol.HasSubstructMatch(thiohydroximate_pattern):
        return True, "Contains valid methylsulfide group in thiohydroximate structure"
    
    # Get molecular weight
    mol_weight = sum([atom.GetMass() for atom in mol.GetAtoms()])
    
    # Filter out very small molecules
    if mol_weight < 60:
        return False, "Molecule too small to be considered a methyl sulfide compound"

    # Verify sulfur connections and environment
    for match in matches:
        s_atom_idx = match[0]  # Index of the sulfur atom
        s_atom = mol.GetAtomWithIdx(s_atom_idx)
        
        # Check degree of sulfur
        if s_atom.GetDegree() != 2:
            return False, "Sulfur has incorrect number of connections"
        
        # Get neighbor atoms
        neighbors = [n for n in s_atom.GetNeighbors() if n.GetIdx() != match[1]]  # Exclude the methyl group
        if not neighbors:
            return False, "Invalid sulfur environment"
            
        neighbor = neighbors[0]
        # If neighbor is aromatic carbon, check if it's part of a heterocycle
        if neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6:
            ring_info = mol.GetRingInfo()
            ring_atoms = ring_info.AtomRings()
            for ring in ring_atoms:
                if neighbor.GetIdx() in ring:
                    # Check if ring contains heteroatoms
                    ring_has_hetero = any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring)
                    if ring_has_hetero:
                        return True, "Contains methylsulfide group attached to heterocyclic ring"
    
    return True, "Contains methylsulfide (S-CH3) group with appropriate structure"