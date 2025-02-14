"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Any fatty acid containing at least one C=C double bond.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is a fatty acid with at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group using SMARTS
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(acid_pattern):
         # Try finding the carboxylic group via atom matching, needed to identify structures where the OH group is represented as O
        acid_atom_match = False
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetBond().GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        for neighbor2 in atom.GetNeighbors():
                            if neighbor2.GetAtomicNum() == 8 and neighbor2.GetBond().GetBondType() == Chem.rdchem.BondType.SINGLE:
                                acid_atom_match = True
                                break
                    if acid_atom_match:
                        break
                if acid_atom_match:
                    break
        if not acid_atom_match:
            return False, "No carboxylic acid group found"

    # Look for at least one C=C double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No C=C double bond found"
    
    # Check for a reasonable number of carbons to be considered a fatty acid (at least 8)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 8:
         return False, "Too few carbons to be a fatty acid."
    
     # Check for a reasonable number of rotatable bonds for a long aliphatic chain
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds < 4:
        return False, "Too few rotatable bonds, not a long fatty acid chain."

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for fatty acid"
    
    # Check number of oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 4:
        return False, "Too many oxygens, not a simple fatty acid."


    return True, "Contains a carboxylic acid group and at least one C=C double bond and has 8 or more carbons, enough rotatable bonds and a reasonable weight."