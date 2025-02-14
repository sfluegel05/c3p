"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:35741 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is any member of the class of chlorobenzenes carrying
    four chloro groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count chlorine atoms
    n_chlorine = sum(atom.GetAtomicNum() == 17 for atom in mol.GetAtoms())
    
    # Tetrachlorobenzene must have exactly 4 chlorine atoms
    if n_chlorine != 4:
        return False, f"Found {n_chlorine} chlorine atoms, need exactly 4"
    
    # Check for benzene ring with 4 chlorine substituents
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1Cl.Cl.Cl.Cl")
    if mol.HasSubstructMatch(benzene_pattern):
        return True, "Contains a benzene ring with 4 chlorine substituents"
    
    # Check for other ring systems with 4 chlorine substituents
    ring_pattern = Chem.MolFromSmarts("*1(*)*2(*)*3(*)*4(*)*5(*)*6*1Cl.Cl.Cl.Cl")
    if mol.HasSubstructMatch(ring_pattern):
        return True, "Contains a ring system with 4 chlorine substituents"
    
    # Check for molecular weight range typical of tetrachlorobenzenes
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 180 and mol_wt < 300:
        return True, "Molecular weight in the range of tetrachlorobenzenes"
    
    return False, "Does not match the expected patterns for tetrachlorobenzenes"