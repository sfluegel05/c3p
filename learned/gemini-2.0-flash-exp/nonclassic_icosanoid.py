"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a non-classic icosanoid based on its SMILES string.
    A non-classic icosanoid is an oxygenated C20 fatty acid that is NOT a leukotriene or prostanoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is a non-classic icosanoid, (False, reason) otherwise
                         or (None, None) if the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
         return None, None

    # Check for 20 Carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 20:
        return False, f"Not a C20 fatty acid: {carbon_count} carbons"

    # Check for at least one oxygen
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms present"

    # Check for leukotriene pattern (conjugated triene, 3 double bonds separated by single bonds)
    # simplified to just look for three double bonds separated by single bonds
    # in a row - more complete substructure search is possible but beyond the scope.
    leukotriene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if mol.HasSubstructMatch(leukotriene_pattern):
         return False, "Likely a Leukotriene (conjugated triene pattern)"

    # Check for prostanoid pattern (cyclopentane ring in the middle)
    prostanoid_pattern = Chem.MolFromSmarts("C1CC[C]1") # simplified cyclopentane
    if mol.HasSubstructMatch(prostanoid_pattern):
        # further verify if cyclopentane is within a C20 chain
        prostanoid_pattern_extended = Chem.MolFromSmarts("C[C]1CC[C]1[C]")
        if mol.HasSubstructMatch(prostanoid_pattern_extended):
            return False, "Likely a Prostanoid (cyclopentane within fatty acid chain)"


    # Check for a carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group"
    
    # Check that it is not just a simple fatty acid
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count == 0:
       hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OH1]"))
       epoxide_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C1OC1"))
       carbonyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C=O"))
       if len(hydroxyl_groups) == 1 and len(epoxide_groups)==0 and len(carbonyl_groups) == 1:
           return False, "Likely just a simple fatty acid."

    return True, "Likely a non-classic icosanoid"