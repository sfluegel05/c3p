"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    This function attempts to identify common characteristics of cannabinoids,
    including the dibenzopyran core for classical cannabinoids, presence of
    oxygen atoms in functional groups or rings, as well as lipid-like chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a tricyclic dibenzopyran ring system (common in classical cannabinoids)
    dibenzopyran_pattern = Chem.MolFromSmarts("c1cc2c3c(cc1)Oc1ccccc1C3CC2")
    has_dibenzopyran = mol.HasSubstructMatch(dibenzopyran_pattern)
    
    # check for the presence of oxygen atoms
    has_oxygen_atoms = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())

    # Check for fatty acid chains.
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    has_fatty_acid_chain = len(fatty_acid_matches) > 0
    
    # Check for an amide bond.
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX2]")
    has_amide_bond = mol.HasSubstructMatch(amide_pattern)

    # Check for ether bond
    ether_pattern = Chem.MolFromSmarts("[OX2]-[CX4]")
    has_ether_bond = mol.HasSubstructMatch(ether_pattern)

    # Check for alcohol functionality
    alcohol_pattern = Chem.MolFromSmarts("[OX2][H]")
    has_alcohol_group = mol.HasSubstructMatch(alcohol_pattern)

     # Check for carbonyl functionality
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    has_carbonyl_group = mol.HasSubstructMatch(carbonyl_pattern)

    # Check for carboxyl functionality
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O[H]")
    has_carboxyl_group = mol.HasSubstructMatch(carboxyl_pattern)

    # Apply the rules
    if has_dibenzopyran and has_oxygen_atoms:
       return True, "Contains dibenzopyran core and oxygen atoms (likely classical cannabinoid)"
    elif has_fatty_acid_chain and (has_amide_bond or has_ether_bond or has_alcohol_group or has_carbonyl_group or has_carboxyl_group):
        return True, "Contains fatty acid chain and oxygen functionality (likely endocannabinoid or related)"
    elif has_oxygen_atoms and has_fatty_acid_chain:
         return True, "Contains oxygen atoms and fatty acid chains"
    else:
       return False, "Does not fit the criteria for a cannabinoid structure."