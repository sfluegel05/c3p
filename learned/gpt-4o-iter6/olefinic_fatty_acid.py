"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    An olefinic fatty acid is characterized by having a carboxylic acid group
    and at least one C=C (carbon-carbon double bond) in an aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group - allow for common variations
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-,H0,H0+]")  # More inclusive of variations
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group or valid variant found"

    # Look for non-cyclic C=C double bonds
    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(cc_double_bond_pattern):
        return False, "No aliphatic carbon-carbon double bond found"
    
    # Check for continuous path connecting a double bond to the carboxylic acid group
    # This helps ensure the double bond is part of the fatty acid-like chain
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("CCC=CCC(=O)O")):
        # General flexible path - require connectivity of chain with constraints
        all_cyes = [bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2]
        for bond in all_cyes:
            if any(atom.GetSymbol() == 'O' for atom in bond.GetBeginAtom().GetNeighbors()):
                break
            if any(atom.GetSymbol() == 'O' for atom in bond.GetEndAtom().GetNeighbors()):
                break
        else:
            return False, "No connection between C=C bonds and carboxylic acid found"

    # Ensure minimum length typically associated with fatty acids
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, "Carbon count too low for typical fatty acids"
    
    # Check excessive functionalization
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in (7, 8, 9, 15, 16))
    if heteroatom_count > 5:
        return False, "Excessive heteroatoms affecting classification consistency"

    return True, "Contains both a carboxylic acid group and at least one aliphatic C=C double bond typical of olefinic fatty acids"