"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33561 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation
    of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplify molecule by adding hydrogens
    mol = Chem.AddHs(mol)
    
    # Define SMARTS patterns
    # Sugar ring (ribose or deoxyribose, furanose ring without specifying stereochemistry)
    sugar_pattern = Chem.MolFromSmarts("[C;!R]=,:[C;!R]1-[O]-[C;!R]-[C;!R]-[C;!R]-1")  # Flexible furanose ring

    # Nitrogenous bases (purine and pyrimidine rings)
    purine_pattern = Chem.MolFromSmarts("c1nc2ncnc-2n1")  # Purine ring
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncncn1")    # Pyrimidine ring

    # Phosphate group (match any phosphate group)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")  # Phosphoric acid

    # Check for sugar ring
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar ring found"

    # Check for nitrogenous base
    base_found = False
    base_patterns = [purine_pattern, pyrimidine_pattern]
    for base_pattern in base_patterns:
        if mol.HasSubstructMatch(base_pattern):
            base_found = True
            break
    if not base_found:
        return False, "No nitrogenous base found"

    # Check for glycosidic bond between sugar and base
    # Find atoms that are connected: sugar anomeric carbon to base nitrogen
    # Define anomeric carbon in sugar (any carbon in sugar ring connected to oxygen)
    anomeric_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 8 and neighbor.IsInRing():
                    anomeric_carbons.append(atom.GetIdx())
    # Check for bond between anomeric carbon and base nitrogen
    glycosidic_bond_found = False
    for anomeric_carbon_idx in anomeric_carbons:
        atom = mol.GetAtomWithIdx(anomeric_carbon_idx)
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 7:  # Nitrogen
                if neighbor.IsInRing():
                    glycosidic_bond_found = True
                    break
        if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
        return False, "No glycosidic bond found between sugar and base"

    # Check for phosphate group attached to sugar
    # Find phosphate groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check if phosphate is connected to sugar at any position
    phosphate_connected_to_sugar = False
    sugar_atom_indices = [atom_idx for match in sugar_matches for atom_idx in match]
    phosphate_atom_indices = [atom_idx for match in phosphate_matches for atom_idx in match]
    for phosphate_idx in phosphate_atom_indices:
        phosphate_atom = mol.GetAtomWithIdx(phosphate_idx)
        neighbors = phosphate_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() in sugar_atom_indices:
                phosphate_connected_to_sugar = True
                break
        if phosphate_connected_to_sugar:
            break
    if not phosphate_connected_to_sugar:
        return False, "Phosphate group not connected to sugar"

    return True, "Molecule is a nucleotide"