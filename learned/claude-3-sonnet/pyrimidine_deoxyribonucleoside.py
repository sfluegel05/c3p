"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base connected to a 
    deoxyribose sugar via an N-glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic 5-membered sugar ring with oxygen
    # More permissive pattern that allows for various substitutions
    sugar_pattern = Chem.MolFromSmarts("[O;R1]1[C;R1][C;R1][C;R1][C;R1]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No furanose sugar ring found"

    # Check for pyrimidine base - look for 6-membered ring with two nitrogens
    # Pattern matches common pyrimidine bases (cytosine, thymine, uracil)
    pyrimidine_pattern = Chem.MolFromSmarts("n1[c;$(C=O)]nc([c;$(C=O)])cc1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # Check for N-glycosidic bond between sugar and base
    # Simplified pattern that just looks for the key connection
    glycosidic_pattern = Chem.MolFromSmarts("[O;R1]1[C;R1][C;R1][C;R1][C;R1]1[N;R1]2[c;R1][n;R1][c;R1][c;R1][c;R1]2")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No N-glycosidic bond between sugar and base"

    # Verify presence of at least one carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing characteristic carbonyl group"

    # Check that the sugar is deoxyribose by looking for correct number of oxygens
    # Get atoms in the sugar ring
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if sugar_matches:
        sugar_atoms = set(sugar_matches[0])
        # Count oxygens connected to sugar ring carbons
        oxygen_count = 0
        for atom_idx in sugar_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # If it's a carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:  # If neighbor is oxygen
                        oxygen_count += 1
        # Deoxyribose should have 3 oxygens (ring O + 2 hydroxyls)
        if oxygen_count > 4:  # Allow one extra for modifications
            return False, "Sugar has too many oxygens to be deoxyribose"
        if oxygen_count < 2:  # Must have at least ring O + 1 hydroxyl
            return False, "Sugar lacks minimum oxygens for deoxyribose"

    # Check total ring count
    rings = mol.GetRingInfo().NumRings()
    if rings > 4:  # Allow up to 4 rings to account for modifications
        return False, "Structure too complex - too many rings"

    return True, "Contains pyrimidine base connected to deoxyribose via N-glycosidic bond"