"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:73333 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define key structural patterns
    patterns = {
        # Deoxyribose sugar with 5' phosphate - more flexible pattern
        # Matches furanose ring connected to phosphate via CH2
        'deoxyribose_5p': '[O,OH1,O-]-[P](=[O])([O,OH1,O-])[O]-[CH2]-[CH]1-[O]-[CH]-[CH]-[CH2]1',
        
        # Common nucleobase patterns
        'pyrimidine': 'c1[c,n]c(=[O,N])[n,c][c,n]1',  # Pyrimidine-based (C, T, U)
        'purine': 'c1nc2[c,n]c[n,c]c2[n,c]1',         # Purine-based (A, G)
        
        # Sugar ring with proper connectivity
        'sugar_ring': '[CH2]1-[CH]-[O]-[CH]-[CH]-[CH2]1'
    }
    
    # Convert patterns to RDKit molecules
    smart_patterns = {}
    for name, pattern in patterns.items():
        smart_mol = Chem.MolFromSmarts(pattern)
        if smart_mol is None:
            return False, f"Invalid SMARTS pattern for {name}"
        smart_patterns[name] = smart_mol
    
    # Check for deoxyribose-phosphate moiety
    if not mol.HasSubstructMatch(smart_patterns['deoxyribose_5p']):
        return False, "No deoxyribose sugar with 5' phosphate found"
    
    # Check for sugar ring with proper connectivity
    if not mol.HasSubstructMatch(smart_patterns['sugar_ring']):
        return False, "No proper deoxyribose ring found"
    
    # Check for nucleobase (either pyrimidine or purine)
    has_pyrimidine = mol.HasSubstructMatch(smart_patterns['pyrimidine'])
    has_purine = mol.HasSubstructMatch(smart_patterns['purine'])
    if not (has_pyrimidine or has_purine):
        return False, "No nucleobase found"
    
    # Count phosphorus atoms - should only have one
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, should be exactly 1"
    
    # Count oxygen atoms (at least 5: 3 from phosphate, 2 from sugar)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Insufficient oxygen atoms for structure"
    
    # Additional validation
    # Verify phosphate is connected to sugar via CH2
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if phos_atoms:
        phos = phos_atoms[0]
        # Get neighbors of phosphorus
        for neighbor in phos.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # oxygen
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetAtomicNum() == 6:  # carbon
                        carbon_neighbors = [n for n in o_neighbor.GetNeighbors()]
                        if len(carbon_neighbors) == 2:  # CH2 group
                            return True, "Contains deoxyribose sugar with 5' phosphate and nucleobase"
    
    return False, "Phosphate not properly connected to 5' position"