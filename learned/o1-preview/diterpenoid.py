"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is any terpenoid derived from a diterpene (C20 skeleton),
    which may be rearranged or modified by the removal of one or more skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carbon atoms (should be at least 15)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Number of carbons ({c_count}) too low for diterpenoid"

    # Check for oxygen atoms (terpenoids typically contain oxygen)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found, terpenoids typically contain oxygen"

    # Get Murcko scaffold to analyze core structure
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return False, "Could not extract molecular scaffold"

    # Count carbons in scaffold
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    if scaffold_c_count < 10:
        return False, f"Scaffold has too few carbons ({scaffold_c_count}) for diterpenoid core"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for diterpenoid"

    # If all checks pass, classify as diterpenoid
    return True, "Molecule matches characteristics of a diterpenoid"