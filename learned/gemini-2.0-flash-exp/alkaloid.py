"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is likely an alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen in a heterocyclic ring
    heterocyclic_nitrogen_pattern = Chem.MolFromSmarts("[#7;R1]~[#6;R1]")
    if not mol.HasSubstructMatch(heterocyclic_nitrogen_pattern):
         return False, "No nitrogen in a heterocyclic ring"


    # Check for exocyclic nitrogen, excluding amides and similar
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.IsInRing():
            # Check if exocyclic nitrogen is part of an amide group
            is_amide_or_similar = False
            for neighbor in atom.GetNeighbors():
                 if neighbor.GetAtomicNum() == 6:
                   for neighbor2 in neighbor.GetNeighbors():
                       if neighbor2.GetAtomicNum() == 8 and neighbor2.GetBondBetweenAtoms(neighbor.GetIdx(),neighbor2.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                            is_amide_or_similar = True
                            break
                 if is_amide_or_similar:
                   break
            if not is_amide_or_similar:
                #Check if nitrogen is only attached to C atoms. If yes, it is considered an exocyclic amine and should not be an alkaloid.
                neighbor_carbon_count = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        neighbor_carbon_count +=1

                if neighbor_carbon_count > 0:
                    all_c_neighbors = True
                    for neighbor in atom.GetNeighbors():
                      if neighbor.GetAtomicNum() !=6:
                        all_c_neighbors = False
                    if all_c_neighbors and neighbor_carbon_count < 3:
                        return False, "Contains exocyclic nitrogen as amine"



    # Check for carboxylic acid groups (-C(=O)O or -C(=O)OH)
    carboxylic_acid_pattern1 = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0]")
    carboxylic_acid_pattern2 = Chem.MolFromSmarts("[CX3](=[OX1])[OH1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern1) or mol.HasSubstructMatch(carboxylic_acid_pattern2):
        return False, "Contains carboxylic acid group"

    #Check for sulfonic acids
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)[OH]")
    if mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "Contains sulfonic acid group"



    return True, "Contains at least one nitrogen in a heterocyclic ring, does not have simple exocyclic amine or carboxylic/sulfonic acid groups"