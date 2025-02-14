"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:38155 dodecanoate ester
A fatty acid ester in which the carboxylic acid component is lauric acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for lauric acid substructure (CCCCCCCCCCCCC(=O)O-)
    lauric_acid = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)[O;X1]")
    lauric_matches = mol.GetSubstructMatches(lauric_acid)

    # Look for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not lauric_matches or not ester_matches:
        return False, "No dodecanoate ester functionality found"

    # Check if ester is connected to lauric acid
    for ester_match in ester_matches:
        ester_oxygen = mol.GetAtomWithIdx(ester_match[1])
        for lauric_atom_idx in lauric_matches[0]:
            lauric_atom = mol.GetAtomWithIdx(lauric_atom_idx)
            if mol.GetBondBetweenAtoms(ester_oxygen.GetIdx(), lauric_atom.GetIdx()):
                # Additional checks
                mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
                if mol_wt < 200 or mol_wt > 500:
                    return False, "Molecular weight outside expected range for dodecanoate esters"

                c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
                if c_count < 12 or o_count < 2:
                    return False, "Insufficient carbon or oxygen atoms for dodecanoate esters"

                return True, "Contains lauric acid ester group"

    return False, "Lauric acid not connected to ester group"