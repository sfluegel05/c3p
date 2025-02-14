"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a flavanone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS with the hydroxyl group at C3
    # The core SMARTS should represent the ring with O-C-C(=O)-C-C-C plus the hydroxyl on the carbon adjacent to the C=O
    flavanone_core_smarts = "[C]1[C](=[O])[C]([OH])[C][C][O]1"
    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core not found"
    
    # Optional checks based on molecular weight, carbon and oxygen count to see if it is close to dihydroflavonols.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for dihydroflavonol"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
      return False, "Too few carbons for dihydroflavonol"
    if o_count < 4:
      return False, "Too few oxygens for dihydroflavonol"


    return True, "Molecule is a dihydroflavonol"