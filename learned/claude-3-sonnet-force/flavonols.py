"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:17773 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is a flavonoid (has benzopyran backbone)
    flavonoid_pattern = Chem.MolFromSmarts("[c1]2[c]([c]3[c]([c]([c]2[o]1)[O])[O])-[c]([c]([c]3=O)[O])")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Not a flavonoid (missing benzopyran backbone)"

    # Check for hydroxy group at position 3 of heterocyclic ring
    position_3_oh_pattern = Chem.MolFromSmarts("[c1]2[c]([c]([c]([o]2)[O])[O])-[c]([c]([c]1=O)[O])[O]")
    if not mol.HasSubstructMatch(position_3_oh_pattern):
        return False, "No hydroxy group at position 3 of heterocyclic ring"

    # Check for additional rings and substituents
    flavonol_pattern = Chem.MolFromSmarts("[c1]2[c]([c]3[c]([c]([c]2[o]1)[O])[O])-[c]([c]([c]3=O)[O])[O]")
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "Missing additional rings or substituents required for flavonols"

    # Check for the correct number of carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 15 or c_count > 25 or o_count < 5 or o_count > 10:
        return False, "Incorrect number of carbons or oxygens for a flavonol"

    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside the typical range for flavonols"

    return True, "Meets structural requirements for a flavonol"