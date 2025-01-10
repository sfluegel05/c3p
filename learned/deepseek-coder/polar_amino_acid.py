"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: CHEBI:26271 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid has an amino group, a carboxyl group, and a side chain capable of forming hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an amino group (-NH2 or -NH3+)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;+0,+1]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Check for the presence of a carboxyl group (-COOH or -COO-)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1H0-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for the presence of a side chain (at least one carbon attached to the alpha carbon)
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H]([NX3])[CX3](=O)")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if not alpha_carbon_matches:
        return False, "No alpha carbon found (not an amino acid)"

    # Check if the side chain is polar (capable of forming hydrogen bonds)
    # Polar groups include -OH, -CONH2, -NH2, -SH, etc.
    polar_pattern = Chem.MolFromSmarts("[OX2H1,OX1H0-,NX3H2,NX3H1,SX2H1]")
    polar_matches = mol.GetSubstructMatches(polar_pattern)
    if len(polar_matches) < 2:  # At least one polar group in the side chain
        return False, "Side chain is not polar (cannot form hydrogen bonds)"

    # Additional check: Ensure the polar group is in the side chain, not part of the backbone
    # The side chain is any atom connected to the alpha carbon that is not part of the amino or carboxyl group
    alpha_carbon_idx = alpha_carbon_matches[0][0]  # Index of the alpha carbon
    side_chain_atoms = set()
    for bond in mol.GetAtomWithIdx(alpha_carbon_idx).GetBonds():
        neighbor = bond.GetOtherAtomIdx(alpha_carbon_idx)
        if not (mol.GetAtomWithIdx(neighbor).GetSymbol() == 'N' or mol.GetAtomWithIdx(neighbor).GetSymbol() == 'C' and mol.GetAtomWithIdx(neighbor).GetDegree() == 3):
            side_chain_atoms.add(neighbor)

    # Check if any polar group is in the side chain
    polar_in_side_chain = False
    for match in polar_matches:
        if match[0] in side_chain_atoms:
            polar_in_side_chain = True
            break

    if not polar_in_side_chain:
        return False, "Polar group is not in the side chain"

    return True, "Contains amino and carboxyl groups with a polar side chain capable of forming hydrogen bonds"