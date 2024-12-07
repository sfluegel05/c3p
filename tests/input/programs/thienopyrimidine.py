from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_thienopyrimidine(smiles: str):
    """
    Determines if a molecule is a thienopyrimidine (pyrimidine ring ortho fused to a 5-membered thiophene ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thienopyrimidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for thienopyrimidine
    thienopyrimidine_pattern = Chem.MolFromSmarts('c1ncnc2scnc12')

    if mol.HasSubstructMatch(thienopyrimidine_pattern):
        return True, "Molecule contains a thienopyrimidine core structure"
    else:
        return False, "Molecule does not contain a thienopyrimidine core structure"

# Example usage:
# print(is_thienopyrimidine('CC(=O)C1=CC=CC(=C1)C2=CSC3=C2C(=NC=N3)SC4=NN=NN4C'))
# print(is_thienopyrimidine('CC1=C(SC2=NC=NC(=C12)NCC(C3=CC=C(C=C3)[N+](=O)[O-])O)C'))

