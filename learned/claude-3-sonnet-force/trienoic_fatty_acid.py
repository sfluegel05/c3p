"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:36195 trienoic fatty acid
A polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count double bonds
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    num_double_bonds = len(double_bonds)
    if num_double_bonds != 3:
        return False, f"Found {num_double_bonds} double bonds, expected 3"

    # Check for long linear carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "No long linear carbon chain found"

    # Check for carboxylic acid group at the end of the chain
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Check if double bonds are in cis or trans configuration
    cis_trans_pattern = Chem.MolFromSmarts("/C=C/")
    cis_trans_matches = mol.GetSubstructMatches(cis_trans_pattern)
    if len(cis_trans_matches) != 3:
        return False, "Double bonds not in expected cis/trans configuration"

    # Check for additional functional groups (optional)
    hydroxy_pattern = Chem.MolFromSmarts("O")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    epoxy_matches = mol.GetSubstructMatches(epoxy_pattern)
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)

    reason = "Contains a long linear carbon chain with three double bonds in cis/trans configuration and a carboxylic acid group"
    if hydroxy_matches:
        reason += ", also contains hydroxy group(s)"
    if epoxy_matches:
        reason += ", also contains epoxy group(s)"
    if hydroperoxy_matches:
        reason += ", also contains hydroperoxy group(s)"

    return True, reason