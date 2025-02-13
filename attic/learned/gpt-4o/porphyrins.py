"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a skeleton of four pyrrole nuclei united through the alpha-positions by four methine groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Porphyrin SMARTS pattern: four pyrrole rings connected by methine bridges
    porphyrin_pattern = Chem.MolFromSmarts("c1nccc1-c2nccc2-c3nccc3-c4nccc4")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "Structure does not match porphyrin macrocycle"

    # Check for potential metal coordination
    # Common porphyrin metals: Fe, Mg, Zn, Co, Ni, etc.
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in ['Fe', 'Mg', 'Zn', 'Co', 'Ni', 'Cu', 'Mn']]
    if len(metal_atoms) == 0:
        return False, "No metal ion commonly associated with porphyrins found"

    return True, "Contains porphyrin macrocycle with potential metal coordination"