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

    # SMARTS patterns for 1,2,4-triazine core
    patterns = [
        # Aromatic form
        "[n]1[n]c[n]cc1",
        # Non-aromatic forms with different bond arrangements
        "N1=NC=NC=C1",
        "N1=NN=CC=C1",
        "N1N=CN=CC1",
        # Fused systems and other variants
        "[n]1[n]c2[n]cc1[n]2",  # For triazolo-triazine
        "[n]1[n]c([n]cc1)*)*",  # For substituted variants
        "[N]1=[N][C](=[N][C]=[C]1)*)*"  # For non-aromatic substituted variants
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
                
                # Count nitrogens in the matched ring
                n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                
                # Verify we have exactly 3 nitrogens
                if n_count != 3:
                    continue
                
                # Get nitrogen positions
                n_positions = []
                for i, atom in enumerate(ring_atoms[:6]):  # Only consider first 6 atoms for base ring
                    if atom.GetAtomicNum() == 7:  # Nitrogen
                        n_positions.append(i)
                
                # Check if the nitrogens form a 1,2,4 pattern
                # Valid patterns are [0,1,3], [1,2,4], [2,3,5], [3,4,0], [4,5,1], [5,0,2]
                valid_patterns = [[0,1,3], [1,2,4], [2,3,5], [3,4,0], [4,5,1], [5,0,2]]
                n_positions = sorted(n_positions)
                
                if n_positions in valid_patterns:
                    # Determine if aromatic
                    is_aromatic = any(atom.GetIsAromatic() for atom in ring_atoms)
                    ring_type = "aromatic" if is_aromatic else "non-aromatic"
                    
                    # Check if part of a fused system
                    is_fused = len(mol.GetRingInfo().AtomRings()) > 1
                    system_type = "fused ring" if is_fused else "single ring"
                    
                    return True, f"Contains {ring_type} 1,2,4-triazine ({system_type}) with nitrogens at positions 1, 2, and 4"

    return False, "No valid 1,2,4-triazine ring found"