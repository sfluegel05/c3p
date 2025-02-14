"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
    The carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific ceramide core pattern - accounts for variations
    # Includes an amide link, two hydroxyls and a glycosidic link to a carbon atom.
    ceramide_pattern = Chem.MolFromSmarts("[CX4][CX4](O)[CX4](O)([NX3][CX3](=[OX1])[#6])[CX4][O]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide/sphingosine core found"
    
    # Check if glycosidic link is at position 1
    # Pattern for glycosidic link at O-1 of sphingosine - [CX4]([OX2][CX4]1)([CX4](O)[CX4](O)N)
    # This pattern looks for the carbon with the glycosidic link attached to the rest of the sphingosine
    glycosidic_pattern = Chem.MolFromSmarts("[CX4]([OX2][CX4]1)([CX4](O)[CX4](O)[NX3])")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage at the correct position"
    
    # Check for long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
       return False, "Chains too short to be sphingolipid"

    # Check molecular weight - glycosphingolipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
         return False, "Molecular weight too low for a glycosphingolipid"

    return True, "Contains a ceramide/sphingosine core with a carbohydrate attached via a glycosidic linkage"