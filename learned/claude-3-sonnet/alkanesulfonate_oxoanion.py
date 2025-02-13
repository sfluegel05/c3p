"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI:35518 alkanesulfonate oxoanion

An alkanesulfonate oxoanion is defined as an alkanesulfonate in which the carbon at position 1 
is attached to R, which can represent hydrogens, a carbon chain, or other groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonate group (S(=O)(=O)[O-]) pattern
    sulfonate_pattern = Chem.MolFromSmarts("[S+](=O)(=O)[O-]")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate group found"
    
    # Look for carbon attached to sulfonate group
    alkanesulfonate_pattern = Chem.MolFromSmarts("[C]-[S+](=O)(=O)[O-]")
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "No carbon attached to sulfonate group"
    
    # Check if the sulfonate group is an oxoanion
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if formal_charge != -1:
        return False, "Formal charge is not -1 (not an oxoanion)"
    
    # Check for additional substituents on the carbon attached to the sulfonate group
    alkanesulfonate_carbon = mol.GetSubstructMatch(alkanesulfonate_pattern)[0]
    alkanesulfonate_carbon_atom = mol.GetAtomWithIdx(alkanesulfonate_carbon)
    substituents = [atom.GetSymbol() for atom in alkanesulfonate_carbon_atom.GetNeighbors()]
    if len(substituents) < 2:
        return False, "No substituents on the carbon attached to the sulfonate group"
    
    return True, "Molecule is an alkanesulfonate oxoanion"