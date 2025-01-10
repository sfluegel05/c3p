"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: organic sulfide compounds (thioethers)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has the structure R-S-R where R is any organic group (not H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check each sulfur atom
    valid_sulfide_count = 0
    for sulfur in sulfur_atoms:
        # Get neighboring atoms
        neighbors = sulfur.GetNeighbors()
        
        # Count number of carbon neighbors
        carbon_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6)
        
        # Check oxidation state (avoid sulfoxides and sulfones)
        explicit_valence = sulfur.GetExplicitValence()
        
        # Check if sulfur has exactly two carbon neighbors and proper valence
        if carbon_neighbors == 2 and explicit_valence == 2:
            valid_sulfide_count += 1
            
        # Check for special cases like thiocyanates (-SCN)
        elif len(neighbors) == 2 and any(n.GetAtomicNum() == 6 for n in neighbors):
            if any(n.GetAtomicNum() == 7 and any(nb.GetAtomicNum() == 6 for nb in n.GetNeighbors()) 
                  for n in neighbors):
                valid_sulfide_count += 1

    if valid_sulfide_count == 0:
        # Check specific patterns that might have been missed
        sulfide_pattern = Chem.MolFromSmarts("[#6]-[#16X2]-[#6]")
        if mol.HasSubstructMatch(sulfide_pattern):
            return True, "Contains R-S-R pattern where R is carbon-containing group"
        return False, "No sulfur atom with two carbon neighbors in correct oxidation state"

    return True, f"Contains {valid_sulfide_count} sulfide group(s) with R-S-R pattern"