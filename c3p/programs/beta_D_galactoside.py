"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside contains a D-galactose ring with beta-configuration at its anomeric center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the D-galactose pyranose ring, without specifying anomeric configuration.
    # This uses absolute stereochemistry markers (@ and @@) to correctly match D-galactose stereochemistry, except at the anomeric carbon
    galactose_ring_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](CO)O[C@@H]1*")
    if not mol.HasSubstructMatch(galactose_ring_pattern):
      return False, "No D-galactose ring found"

    # Define a SMARTS pattern for the beta-D-galactoside, now specifying the beta configuration at the anomeric carbon using [C@H]
    beta_galactose_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@@H](O)[C@@H](O)[C@H](CO)O[C@@H]1*")
    if mol.HasSubstructMatch(beta_galactose_pattern):
        return True, "Contains a beta-D-galactoside moiety"

    return False, "Not a beta-D-galactoside, based on anomeric carbon stereochemistry"