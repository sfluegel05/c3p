"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for monosaccharide units and glycosidic bonds
    monosaccharide_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O')
    glycosidic_bond_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O[C@@H]2[C@@H](CO)O[C@@H]([C@H](O)[C@H]2O)O)[C@@H](O)[C@H](O)[C@@H]1')
    
    # Count monosaccharide units
    num_monosaccharides = len(mol.GetSubstructMatches(monosaccharide_pattern))
    
    # Count glycosidic bonds
    num_glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_bond_pattern))
    
    # Check for continuous chain of monosaccharide units
    continuous_chain = check_continuous_chain(mol, glycosidic_bond_pattern)
    
    # Check criteria for tetrasaccharide
    if num_monosaccharides == 4 and num_glycosidic_bonds >= 3 and continuous_chain:
        return True, "Contains 4 monosaccharide units connected by at least 3 glycosidic bonds"
    else:
        return False, f"Found {num_monosaccharides} monosaccharide units and {num_glycosidic_bonds} glycosidic bonds, but the structure does not match a tetrasaccharide"

def check_continuous_chain(mol, glycosidic_bond_pattern):
    """
    Checks if the monosaccharide units form a continuous chain connected by glycosidic bonds.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
        glycosidic_bond_pattern (rdkit.Chem.rdchem.Mol): SMARTS pattern for glycosidic bonds

    Returns:
        bool: True if the monosaccharide units form a continuous chain, False otherwise
    """
    atom_indices = set()
    for match in mol.GetSubstructMatches(glycosidic_bond_pattern):
        for idx in match:
            atom_indices.add(idx)
    
    num_atoms = mol.GetNumAtoms()
    if len(atom_indices) == num_atoms:
        return True
    else:
        return False

# Example usage
print(is_tetrasaccharide('OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O[C@@H]4[C@@H](CO)OC(O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O'))
# Output: (True, 'Contains 4 monosaccharide units connected by at least 3 glycosidic bonds')