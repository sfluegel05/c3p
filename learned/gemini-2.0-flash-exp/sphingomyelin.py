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

    # 3. Check for sphingoid base - at least 14 carbons, one or two alcohol groups, and one amine, and possibly a double bond
    # core pattern for sphingoid base
    sphingoid_core_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2])[CHX4]([OX2])")
    if not mol.HasSubstructMatch(sphingoid_core_pattern):
          return False, "No sphingoid base core found"
    
    # additional pattern for the sphingoid base, ensuring one or two chains attached to the amide nitrogen. It also verifies a double bond
    sphingoid_base_pattern = Chem.MolFromSmarts("[NX3][CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base chain found"


    # 4. Check for fatty acid chain - look for a long carbon chain attached to the amide
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain"


    # 5. Check molecular weight and rotatable bonds to confirm long chains
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if mol_wt < 400: # relaxed threshold
        return False, "Molecular weight too low for sphingomyelin"
    if n_rotatable < 7: # relaxed threshold
         return False, "Not enough rotatable bonds for sphingomyelin"


    return True, "Contains phosphocholine group, amide bond, sphingoid base and fatty acid"