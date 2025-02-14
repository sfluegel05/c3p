"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is likely an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an endocannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"


    # 1. Check for long carbon chain with double bonds
    # fatty_acid_chain pattern: long carbon chain with a few double bonds
    fatty_acid_pattern = Chem.MolFromSmarts("C~C=C~C~C=C~C") 
    if not mol.HasSubstructMatch(fatty_acid_pattern):
       fatty_acid_pattern = Chem.MolFromSmarts("C~C=C~C~C=C~C~C=C~C")
       if not mol.HasSubstructMatch(fatty_acid_pattern):
           fatty_acid_pattern = Chem.MolFromSmarts("C~C=C~C~C=C~C~C=C~C~C=C~C")
           if not mol.HasSubstructMatch(fatty_acid_pattern):
                fatty_acid_pattern = Chem.MolFromSmarts("C~C~C~C~C~C~C~C~C~C") #also allow for very long saturated chains
                if not mol.HasSubstructMatch(fatty_acid_pattern):
                     return False, "No characteristic fatty acid chain found"


    # 2. Check for polar head group (ethanolamide or glycerol or similar) and linkage (ester or amide)
    ethanolamide_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    # simplified ester and amide for link to long chain
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ether_pattern = Chem.MolFromSmarts("COC")

    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)



    if not (has_ethanolamide or has_glycerol) and not (has_ester or has_amide or has_ether) :
        return False, "No polar head group and linkage found"



    #check for connection
    if has_amide:
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        for match in amide_matches:
            carbon_atoms_linked_to_amide_C = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
            if len(carbon_atoms_linked_to_amide_C) > 0:
                break
        else:
            return False, "No long chain attached to amide"
    if has_ester:
       ester_matches = mol.GetSubstructMatches(ester_pattern)
       for match in ester_matches:
          carbon_atoms_linked_to_ester_C = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
          if len(carbon_atoms_linked_to_ester_C) > 0:
            break
       else:
            return False, "No long chain attached to ester"
    if has_ether:
        ether_matches = mol.GetSubstructMatches(ether_pattern)
        for match in ether_matches:
             carbon_atoms_linked_to_ether_O = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in mol.GetAtomWithIdx(match[1]).GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
             if len(carbon_atoms_linked_to_ether_O) > 0:
                 break
        else:
           return False, "No long chain attached to ether"
    # If all conditions met, it is likely an endocannabinoid
    return True, "Likely an endocannabinoid based on long chain, polar head group, and linker."