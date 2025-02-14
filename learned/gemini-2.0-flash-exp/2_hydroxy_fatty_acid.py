"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group
    on the carbon alpha to the carboxyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(acid_pattern):
         return False, "No carboxylic acid group found"

    # Check for 2-hydroxy group (hydroxy on the alpha carbon)
    # Relaxed pattern to allow for any carbon at the alpha position
    hydroxy_pattern = Chem.MolFromSmarts("[CX4,CX3]([OH])-[CX3](=[O])[O,OH]") 
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxyl group at the alpha position found"
    

    # Check for number of rotatable bonds for chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4: # Fatty acid chains are typically long
      return False, "Chain too short or not flexible enough to be a fatty acid"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecule too small to be a fatty acid"

    # Check for ring system - fatty acids are usually not cyclic.
    ring_info = mol.GetRingInfo()
    if ring_info.IsAtomInRing(mol.GetAtoms()[0].GetIdx()):
        return False, "Molecule is cyclic"


    return True, "Contains a carboxylic acid group with a hydroxyl group on the alpha carbon and a carbon chain"