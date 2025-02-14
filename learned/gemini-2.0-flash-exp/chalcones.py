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
    chalcone_core_pattern = Chem.MolFromSmarts("[c;R]~[C]=[C]~[C](=[O])~[c;R]")
    core_matches = mol.GetSubstructMatches(chalcone_core_pattern)

    if not core_matches:
        return False, "Core chalcone structure (Ar-CH=CH-C(=O)-Ar) not found"


    # Check for aromatic rings attached to the core
    for match in core_matches:
        atom_indices = list(match)
        aryl1_index = atom_indices[0]
        aryl2_index = atom_indices[4]
        aryl1_atom = mol.GetAtomWithIdx(aryl1_index)
        aryl2_atom = mol.GetAtomWithIdx(aryl2_index)
        
        # Check if the terminal groups are aromatic
        if not (aryl1_atom.GetIsAromatic() and aryl2_atom.GetIsAromatic()):
          return False, "At least one of the terminal groups aren't aromatic"


        # Check if the aromatic rings are directly attached
        
        neighbor1 = mol.GetAtomWithIdx(atom_indices[1])
        neighbor2 = mol.GetAtomWithIdx(atom_indices[3])
        
        bond1 = mol.GetBondBetweenAtoms(aryl1_index, neighbor1.GetIdx())
        bond2 = mol.GetBondBetweenAtoms(aryl2_index, neighbor2.GetIdx())
        
        if not (bond1 and bond1.GetBondType() == Chem.BondType.SINGLE):
            return False, "Aryl group 1 is not directly attached to the core"
        
        if not (bond2 and bond2.GetBondType() == Chem.BondType.SINGLE):
            return False, "Aryl group 2 is not directly attached to the core"
        

    # Additional substitution check - OH, OMe or prenyl (not enforced)
    substituent_pattern = Chem.MolFromSmarts("[OH,Oc,CC=C(C)C]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)


    # Check molecular weight - chalcones usually > 200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a chalcone"

    return True, "Contains the core chalcone structure (Ar-CH=CH-C(=O)-Ar) and two aromatic rings"