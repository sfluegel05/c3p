"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a generalized quinone core with at least 2 double bonds and 2 =O groups (allowing for substitutions)
    quinone_pattern = Chem.MolFromSmarts("[C;$(C=O)]1=[C,c][C,c](=[C,c][C,c](=[C,c]1)[C;$(C=O)])=[C,c]") #General six membered ring with 2 carbonyls, and 2 double bonds
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone core found"

    # Define a pattern for a prenyl chain attached to the quinone core
    #The chain can have at least two isoprene units connected by a single bond
    # We'll use multiple SMARTS patterns to capture different variations
    prenyl_patterns = [
      Chem.MolFromSmarts("[C;$(C=O)]~[C,c]!@[CH2]C(=C)C~[CH2]C=C[CX4]"), #Quinone core - C - C=C -C -C=C-
      Chem.MolFromSmarts("[C;$(C=O)]~[C,c]!@[CH2]C(=C)C~[CH2]C(=C)C~[CH2]C=C[CX4]"),#Quinone core - C - C=C -C -C=C -C-C=C-
      Chem.MolFromSmarts("[C;$(C=O)]~[C,c]!@[CH2]C(=C)C~[CH2]C(=C)C~[CH2]C(=C)C~[CH2]C=C[CX4]"), #Quinone core - C - C=C -C -C=C -C -C=C -C-C=C-
      Chem.MolFromSmarts("[C;$(C=O)]~[C,c]!@[CH2]C(=C)C~[CH2]C=C[CX4]~[CH2]C(=C)C~[CH2]C=C[CX4]"), #Quinone core - C - C=C -C -C=C-C-C=C -C-C=C-
      Chem.MolFromSmarts("[C;$(C=O)]~[C,c]!@[CH2]C(=C)C~[CH2]C(=C)C~[CH2]C(=C)C~[CH2]C(=C)C~[CH2]C=C[CX4]")#Quinone core - C - C=C -C -C=C -C -C=C -C -C=C-C-C=C-

      ]
      
    found_prenyl = False
    for pattern in prenyl_patterns:
        if mol.HasSubstructMatch(pattern):
            found_prenyl = True
            break #Exit loop as soon as at least one prenyl chain is found

    if not found_prenyl:
          return False, "No polyprenyl side chain found attached to the quinone"

    
    #Check for minimum number of carbons in the molecule to filter out false positives
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 13: #Arbitrary number based on the examples. A molecule should have at least 13 carbons to be a prenylquinone
        return False, "Too few carbons for prenylquinone"

    return True, "Contains a quinone core and a polyprenyl side chain"