"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
    This term is commonly used only for those containing more than ten monosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more accurate pattern for a monosaccharide unit
    monosaccharide_pattern = Chem.MolFromSmarts("[C;H1][C;H1][C;H1][C;H1][C;H1][O;H1]")
    
    # Find all matches of the monosaccharide pattern
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    # Count the number of distinct monosaccharide units
    num_monosaccharides = len(set(match[0] for match in monosaccharide_matches))
    
    # Check for glycosidic bonds
    glycosidic_bond_pattern = Chem.MolFromSmarts("[O;H0][C;H1][C;H1][O;H0]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check if the number of monosaccharide units is greater than 10
    # and there are glycosidic bonds
    # and molecular weight is high
    if (num_monosaccharides > 10 and 
        len(glycosidic_bond_matches) > 0 and 
        mol_wt > 5000):
        return True, f"Contains {num_monosaccharides} monosaccharide units with glycosidic bonds and high molecular weight ({mol_wt:.1f} Da)"
    else:
        return False, f"Contains {num_monosaccharides} monosaccharide units, which is not more than 10, or lacks glycosidic bonds, or has low molecular weight ({mol_wt:.1f} Da)"