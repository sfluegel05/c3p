"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is defined as a 1,3-diphenylpropenone (benzylideneacetophenone) and its
    derivatives formed by substitution, Ar-CH=CH-C(=O)-Ar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core chalcone SMARTS pattern
    chalcone_core_pattern = Chem.MolFromSmarts("[c]~[C]=[C]~[C](=[O])~[c]")
    core_matches = mol.GetSubstructMatches(chalcone_core_pattern)

    if not core_matches:
        return False, "Core chalcone structure (Ar-CH=CH-C(=O)-Ar) not found"

    # Check for 2 different aryl groups attached to the core
    for match in core_matches:
      atom_indices = list(match)
      aryl1_index = atom_indices[0]
      aryl2_index = atom_indices[4]
      aryl1_atom = mol.GetAtomWithIdx(aryl1_index)
      aryl2_atom = mol.GetAtomWithIdx(aryl2_index)
      
      aryl1_ring = aryl1_atom.IsInRing()
      aryl2_ring = aryl2_atom.IsInRing()
      
      if not (aryl1_ring and aryl2_ring):
        return False, "At least one of the terminal groups aren't aromatic"
        
      aryl1_smiles = Chem.MolFragmentToSmiles(mol, [x.GetIdx() for x in mol.GetAtomNeighbors(aryl1_atom)])
      aryl2_smiles = Chem.MolFragmentToSmiles(mol, [x.GetIdx() for x in mol.GetAtomNeighbors(aryl2_atom)])
      if aryl1_smiles == aryl2_smiles:
           return False, "Terminal aryl groups must be different"

    # Additional substitution check - OH, OMe or prenyl
    substituent_pattern = Chem.MolFromSmarts("[OH,Oc,CC=C(C)C]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)

    # Check molecular weight - chalcones usually > 200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a chalcone"

    return True, "Contains the core chalcone structure (Ar-CH=CH-C(=O)-Ar) and two aromatic rings"