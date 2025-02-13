"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35813 carboxamidine

Carboxamidines are compounds having the structure RC(=NR)NR2. The term is used
as a suffix in systematic nomenclature to denote the -C(=NH)NH2 group including
its carbon atom.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxamidine substructure
    carboxamidine_pattern = Chem.MolFromSmarts("[NX3][CX3]([NX3])=[NX2]")
    matches = mol.GetSubstructMatches(carboxamidine_pattern)

    if not matches:
        return False, "No carboxamidine substructure found"

    # Check if the matched substructure is indeed a carboxamidine group
    for match in matches:
        carbon_idx = match[1]
        nitrogen_idxs = [match[0], match[2]]

        # Check if the carbon atom has a double bond to one nitrogen
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if sum(bond.GetBondTypeAsDouble() == 2 for bond in carbon_atom.GetBonds()) != 1:
            continue

        # Check if the other two nitrogen atoms are connected by a single bond
        nitrogen_atoms = [mol.GetAtomWithIdx(idx) for idx in nitrogen_idxs]
        if any(bond.GetBondTypeAsDouble() != 1 for atom1 in nitrogen_atoms for atom2 in nitrogen_atoms for bond in atom1.GetBonds(atom2)):
            continue

        # Check if the carbon atom has at least one substituent
        if sum(bond.GetBondTypeAsDouble() == 1 for bond in carbon_atom.GetBonds()) < 1:
            return False, "Carboxamidine group must have at least one substituent"

        # Check if the two nitrogen atoms have at least one substituent each
        for atom in nitrogen_atoms:
            if sum(bond.GetBondTypeAsDouble() == 1 for bond in atom.GetBonds()) < 1:
                return False, "Carboxamidine group requires at least one substituent on each nitrogen atom"

        # If all checks pass, it is a valid carboxamidine group
        return True, "Molecule contains the carboxamidine substructure R-C(=NR)-NR2"

    # If no match passes all checks, it is not a carboxamidine
    return False, "No valid carboxamidine substructure found"