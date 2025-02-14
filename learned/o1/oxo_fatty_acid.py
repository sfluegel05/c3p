"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Assume fatty acids have at least one carboxylic acid group
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)

    # Identify aldehyde groups (-CHO)
    aldehyde = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde)
    
    # Identify ketone groups (>C=O)
    ketone = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone)
    
    # Exclude carboxylic acid carbon from aldehyde and ketone matches
    oxo_matches = []
    for match in aldehyde_matches + ketone_matches:
        if not any(idx in carboxylic_acid_atoms for idx in match):
            oxo_matches.append(match)
    
    if not oxo_matches:
        return False, "No additional aldehydic or ketonic group found"

    # Identify the carbon chains connected to carboxylic acid groups
    # and check if oxo groups are part of these chains
    for carboxylic_acid_match in carboxylic_acid_matches:
        carboxylic_carbon_idx = carboxylic_acid_match[0]
        # Perform a BFS to find the longest carbon chain
        visited = set()
        queue = [(carboxylic_carbon_idx, [])]
        chain_atoms = set()
        while queue:
            atom_idx, path = queue.pop(0)
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            # Consider carbon atoms and allow for double bonds and rings
            if atom.GetAtomicNum() == 6:
                chain_atoms.add(atom_idx)
                neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                for nbr_idx in neighbors:
                    if nbr_idx not in visited:
                        queue.append((nbr_idx, path + [atom_idx]))
            elif atom_idx == carboxylic_carbon_idx:
                # Include the immediate neighbor (could be oxygen in carboxyl group)
                neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                for nbr_idx in neighbors:
                    if nbr_idx not in visited:
                        queue.append((nbr_idx, path + [atom_idx]))
        # Check chain length (minimum 4 carbons)
        if len(chain_atoms) < 4:
            continue  # Try next carboxylic acid group
        # Check if any oxo group is part of the chain
        oxo_on_chain = False
        for match in oxo_matches:
            if any(idx in chain_atoms for idx in match):
                oxo_on_chain = True
                break
        if oxo_on_chain:
            # Additional check for unwanted functional groups
            unwanted_groups = [
                Chem.MolFromSmarts("N"),  # Exclude amines/amides
                Chem.MolFromSmarts("C(=O)N"),  # Exclude amides
                Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Exclude esters
                Chem.MolFromSmarts("C(=O)O[C,H]"),  # Exclude esters
                Chem.MolFromSmarts("c"),  # Exclude aromatic rings
                Chem.MolFromSmarts("[#6][OH]"),  # Exclude alcohols (should allow one)
            ]
            found_unwanted = False
            for group in unwanted_groups:
                if mol.HasSubstructMatch(group):
                    found_unwanted = True
                    break
            if found_unwanted:
                continue  # Skip this molecule
            else:
                return True, "Valid oxo fatty acid: contains carboxylic acid and additional oxo group on the fatty acid chain"
    return False, "Does not meet criteria for an oxo fatty acid"