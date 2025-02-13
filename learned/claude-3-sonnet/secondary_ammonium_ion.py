"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:51953 secondary ammonium ion
An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for positively charged atoms
    charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0]
    if len(charged_atoms) == 0:
        return False, "No charged atoms found"
    elif len(charged_atoms) > 1:
        return False, "Multiple charged atoms found"

    # Get the charged atom
    charged_atom = charged_atoms[0]

    # Check if the charged atom is nitrogen with positive charge and sp3 hybridization
    if charged_atom.GetAtomicNum() != 7 or charged_atom.GetFormalCharge() != 1 or charged_atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False, "Charged atom is not a protonated secondary amine"

    # Check if the nitrogen has exactly two substituents
    if len(charged_atom.GetNeighbors()) != 3:
        return False, "Charged nitrogen does not have exactly two substituents"

    # Check if the molecule has a net charge of +1
    mol_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if mol_charge != 1:
        return False, "Molecule does not have a net charge of +1"

    return True, "Molecule contains a protonated secondary amine"