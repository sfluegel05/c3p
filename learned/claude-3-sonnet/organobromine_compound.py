"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:33257 organobromine compound
A compound containing at least one carbon-bromine bond.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon-bromine bonds
    carbon_bromine_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("[Br;X3]~[#6]"))
    if not carbon_bromine_bonds:
        return False, "No carbon-bromine bonds found"

    # Check for aromaticity
    aromatic_rings = mol.GetAromaticRings()
    for ring in aromatic_rings:
        for bond in ring:
            atom1 = mol.GetAtomWithIdx(bond[0])
            atom2 = mol.GetAtomWithIdx(bond[1])
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 35 and bond.GetBondType() == Chem.BondType.AROMATIC:
                return True, "Contains an aromatic carbon-bromine bond"

    # Check bond orders
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 35 and bond.GetBondType() == Chem.BondType.SINGLE:
            return True, "Contains a single carbon-bromine bond"

    return False, "No valid carbon-bromine bonds found"