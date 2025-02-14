"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelins are phospholipids with a sphingoid base, a fatty acid, and a phosphocholine headgroup.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[OX2]-[P](=[OX1])-[OX2]-[CX4]-[CX4]-[N+]([CX4])([CX4])[CX4]")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # 2. Check for Amide linkage
    amide_pattern = Chem.MolFromSmarts("[NX3]-[CX3](=[OX1])")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # 3. Check for sphingoid base - improved pattern - two alcohols and an amine connected to amide and long chains
    # Modified core pattern, focusing on the three carbon chain with two alcohols and a double bond.
    # Modified, includes a third alcohol
    sphingoid_core_pattern = Chem.MolFromSmarts("[NX3][CX4]([OX2])[CX4]([OX2])([CX4])[CX4,CX3]=[CX4,CX3]") #this pattern contains double bond
    sphingoid_core_pattern_2 = Chem.MolFromSmarts("[NX3][CX4]([OX2])[CX4]([OX2])([OX2])[CX4]") #this pattern includes the 3 alcohol, when there is a 4-OH

    if not mol.HasSubstructMatch(sphingoid_core_pattern) and not mol.HasSubstructMatch(sphingoid_core_pattern_2):
          return False, "No sphingoid base core found"
    
    #ensure that the nitrogen is directly attached to a carbonyl group
    sphingoid_amide_pattern = Chem.MolFromSmarts("[NX3]([CX4,CX3])-[CX3](=[OX1])")
    if not mol.HasSubstructMatch(sphingoid_amide_pattern):
        return False, "Sphingoid base not connected to a carbonyl"

    # 4. Check for fatty acid chain -  allow [CX4] or [CX3] next to carbonyl, and ensure it is connected to the sphingoid amine
    fatty_acid_pattern = Chem.MolFromSmarts("[NX3]-[CX3,CX4](=[OX1])-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain found"

    # 5. Check molecular weight and rotatable bonds to confirm long chains
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if mol_wt < 500: # relaxed threshold
        return False, "Molecular weight too low for sphingomyelin"
    if n_rotatable < 8: # relaxed threshold
         return False, "Not enough rotatable bonds for sphingomyelin"

    return True, "Contains phosphocholine group, amide bond, sphingoid base and fatty acid"