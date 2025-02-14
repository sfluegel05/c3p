"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:37585 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for common nucleobase analogue structural features
    purine_pattern = Chem.MolFromSmarts("c1nc2[nH]cnc2[nH]c1")
    pyrimidine_pattern = Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1")
    amino_pattern = Chem.MolFromSmarts("[NH2,NH]")
    oxo_pattern = Chem.MolFromSmarts("C=O")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    alkyl_pattern = Chem.MolFromSmarts("[CX4H3]")
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    
    # Check for purine or pyrimidine ring
    if not mol.HasSubstructMatch(purine_pattern) and not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No purine or pyrimidine ring found"
    
    # Check for amino, oxo, hydroxyl, alkyl, and halogen functional groups
    has_amino = mol.HasSubstructMatch(amino_pattern)
    has_oxo = mol.HasSubstructMatch(oxo_pattern)
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_alkyl = mol.HasSubstructMatch(alkyl_pattern)
    has_halogen = mol.HasSubstructMatch(halogen_pattern)
    if not (has_amino or has_oxo or has_hydroxyl or has_alkyl or has_halogen):
        return False, "Missing typical functional groups found in nucleobase analogues"
    
    # Check for absence of carboxylic acids and esters
    if mol.HasSubstructMatch(carboxyl_pattern) or mol.HasSubstructMatch(ester_pattern):
        return False, "Molecule contains carboxylic acid or ester groups, which are uncommon in nucleobase analogues"
    
    # Check for typical size range of nucleobase analogues
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    if n_heavy_atoms < 5 or n_heavy_atoms > 25:
        return False, "Size deviates too much from typical nucleobases"
    
    # Check for aromaticity and planarity
    aromatic_rings = Chem.GetAromaticRings(mol)
    if not aromatic_rings or not AllChem.PlanarityAtomsSelfIter(mol):
        return False, "Molecule is not sufficiently aromatic or planar"
    
    # Check molecular weight
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 50 or mol_wt > 300:
        return False, "Molecular weight deviates too much from typical nucleobase analogues"
    
    # If all conditions are met, classify as nucleobase analogue
    return True, "Contains purine or pyrimidine ring and typical functional groups of nucleobase analogues"