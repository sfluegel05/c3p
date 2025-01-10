"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:38047 1,2,4-triazine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine ring structure.
    1,2,4-triazine has nitrogen atoms at positions 1, 2, and 4 of a 6-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 1,2,4-triazine ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to standardize aromatic bonds
    Chem.Kekulize(mol, clearAromaticFlags=True)
    
    # SMARTS patterns for 1,2,4-triazine core
    # More specific patterns that explicitly define the 1,2,4-triazine structure
    patterns = [
        # Basic 1,2,4-triazine pattern
        "[N]1=[N][C](=[N][C]=[C]1)",
        # Aromatic form
        "n1nccnc1",
        # Common tautomeric forms
        "[NH]1[N]=[C]([N]=[C][C]1)",
        "[N]1[NH][C](=[N][C]=[C]1)",
        # Reduced forms
        "[NH]1[NH][C](=[N][C]=[C]1)",
    ]

    for pattern in patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if substructure is None:
            continue
            
        if mol.HasSubstructMatch(substructure):
            matches = mol.GetSubstructMatches(substructure)
            for match in matches:
                # Get the matched atoms
                ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
                
                # Verify ring size
                if len(ring_atoms) != 6:
                    continue
                
                # Count nitrogens and verify their positions
                n_positions = []
                for i, atom in enumerate(ring_atoms):
                    if atom.GetAtomicNum() == 7:  # Nitrogen
                        n_positions.append(i)
                
                # Must have exactly 3 nitrogens
                if len(n_positions) != 3:
                    continue
                
                # Check nitrogen positions
                # In a 1,2,4-triazine, nitrogens must be at positions 1,2,4 or equivalent
                n_positions = sorted(n_positions)
                valid_patterns = [
                    [0,1,3],  # Standard 1,2,4 pattern
                    [1,2,4],  # Rotated by 1
                    [2,3,5],  # Rotated by 2
                    [3,4,0],  # Rotated by 3
                    [4,5,1],  # Rotated by 4
                    [5,0,2]   # Rotated by 5
                ]
                
                if n_positions not in valid_patterns:
                    continue
                
                # Verify ring planarity and conjugation
                is_planar = True
                for atom in ring_atoms:
                    # Check hybridization
                    if atom.GetHybridization() not in [Chem.HybridizationType.SP2, 
                                                     Chem.HybridizationType.AROMATIC]:
                        is_planar = False
                        break
                
                if not is_planar:
                    continue
                
                # Determine if part of a fused system
                ring_info = mol.GetRingInfo()
                atom_rings = ring_info.AtomRings()
                is_fused = False
                ring_atoms_set = set(match)
                for ring in atom_rings:
                    if ring_atoms_set != set(ring) and ring_atoms_set.intersection(set(ring)):
                        is_fused = True
                        break
                
                ring_type = "fused" if is_fused else "single"
                return True, f"Contains 1,2,4-triazine {ring_type} ring with correct nitrogen positions"

    return False, "No valid 1,2,4-triazine ring found"