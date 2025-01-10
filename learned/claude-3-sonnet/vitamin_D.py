"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize scoring system
    score = 0
    features = []

    # Characteristic vitamin D triene patterns (more specific)
    triene_patterns = [
        Chem.MolFromSmarts("C=CC=CC=C"),  # Basic triene
        Chem.MolFromSmarts("[CH2][CH2]C=CC=CC=C"),  # With specific context
        Chem.MolFromSmarts("C1CCC(=CC=CC=C)CC1")  # Ring-connected triene
    ]
    
    triene_found = any(mol.HasSubstructMatch(pattern) for pattern in triene_patterns)
    if triene_found:
        score += 2
        features.append("triene system")

    # Seco-steroid skeleton patterns
    seco_patterns = [
        Chem.MolFromSmarts("C1CC[C@H]2[C@@H]1CC[C@H]2C"),  # CD ring system
        Chem.MolFromSmarts("C1CC(=C)CC1"),  # A-ring pattern
        Chem.MolFromSmarts("C1CCC2(C)C1CCC2")  # Modified CD rings
    ]
    
    seco_matches = sum(1 for pattern in seco_patterns if mol.HasSubstructMatch(pattern))
    if seco_matches >= 1:
        score += seco_matches
        features.append("seco-steroid skeleton")

    # Hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches >= 1:
        score += min(hydroxyl_matches, 3)  # Cap at 3
        features.append(f"{hydroxyl_matches} hydroxyl groups")

    # Side chain patterns
    side_chain_patterns = [
        Chem.MolFromSmarts("CCCC(C)C"),  # Standard
        Chem.MolFromSmarts("CCCC(O)(C)C"),  # 25-hydroxy
        Chem.MolFromSmarts("CCC(O)C(C)C")  # Modified
    ]
    
    if any(mol.HasSubstructMatch(pattern) for pattern in side_chain_patterns):
        score += 1
        features.append("characteristic side chain")

    # Basic structural requirements
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if 20 <= c_count <= 35:
        score += 1
        features.append(f"{c_count} carbons")

    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if 2 <= ring_count <= 4:
        score += 1
        features.append(f"{ring_count} rings")

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 340 <= mol_wt <= 700:  # Expanded range
        score += 1
        features.append(f"MW {mol_wt:.1f}")

    # Require minimum score for classification
    if score < 6:
        return False, f"Insufficient vitamin D characteristics (score {score})"

    # Check for common non-vitamin D structures
    non_vit_d_patterns = [
        Chem.MolFromSmarts("C1=CC=CC=C1C1=CC=CC=C1"),  # Biphenyl
        Chem.MolFromSmarts("C1=CC=C2C(=O)C=CC2=C1"),   # Naphthoquinone
        Chem.MolFromSmarts("C=CC=CC=CC=CC=CC=CC=C")    # Extended polyene
    ]
    
    if any(mol.HasSubstructMatch(pattern) for pattern in non_vit_d_patterns):
        score -= 2

    if score < 6:
        return False, "Structure more characteristic of non-vitamin D compound"

    return True, "Contains vitamin D characteristics: " + ", ".join(features)