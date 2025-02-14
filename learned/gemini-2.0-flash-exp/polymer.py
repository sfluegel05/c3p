"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re

def is_polymer(smiles: str):
    """
    Determines if a molecule is likely a polymer based on its SMILES string.
    This is a heuristic approach and may not be 100% accurate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a polymer, False otherwise
        str: Reason for classification
    """

    if smiles is None or not isinstance(smiles,str) or len(smiles) == 0:
        return False, "Invalid SMILES string"

    # Handle mixtures and salts
    smiles_parts = smiles.split('.')
    if len(smiles_parts) > 1:
        is_polymer_mixture = False
        reasons = []
        for part in smiles_parts:
             mol = Chem.MolFromSmiles(part)
             if mol is None:
                 reasons.append("Invalid SMILES string")
                 continue
             is_polymer_part, reason = _check_polymer(mol)
             if is_polymer_part:
                is_polymer_mixture = True
                reasons.append(reason)
        if is_polymer_mixture:
             return True, "Mixture containing polymer: " + ", ".join(reasons)

        return False, "Mixture of small molecules, not a polymer"
    
    # If not a mixture, check as single molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    return _check_polymer(mol)

def _check_polymer(mol):
    """
    Helper function to check if a single molecule is a polymer
    """

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Check for long chains and high MW
    if mol_wt > 1200 and num_rotatable_bonds > 25:
        return True, "High molecular weight and long chain, likely a polymer"

    # Check for repeating patterns
    chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]") # 8 carbon chain
    if mol.HasSubstructMatch(chain_pattern) and mol_wt > 600:
        return True, "Contains long carbon chain, likely a polymer"

    # Check for repeating isoprene units
    isoprene_pattern = Chem.MolFromSmarts("CC(C)=CC")
    if mol.HasSubstructMatch(isoprene_pattern) and mol_wt > 300 and num_rotatable_bonds > 10 :
            return True, "Contains repeating isoprene units, likely a polymer"
    
    # Check for repeating units using a more generic pattern
    repeating_unit_pattern = Chem.MolFromSmarts("[*]~[*]~[*]") #check for at least 3 repeating units
    num_repeating_units = len(mol.GetSubstructMatches(repeating_unit_pattern))

    if num_repeating_units > 5 and mol_wt > 400:
       return True, "Contains repeating units, potentially a polymer"


    # Check for phosphate, sulfate, carboxylate and other potential functional groups (biopolymers)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)(O)")
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    num_phosphates = len(mol.GetSubstructMatches(phosphate_pattern))
    num_sulfates = len(mol.GetSubstructMatches(sulfate_pattern))
    num_carboxylates = len(mol.GetSubstructMatches(carboxylate_pattern))
    num_amides = len(mol.GetSubstructMatches(amide_pattern))
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    num_esters = len(mol.GetSubstructMatches(ester_pattern))

    if (num_phosphates + num_sulfates + num_carboxylates) > 2 and mol_wt > 400: #biopolymers
        return True, "Contains multiple phosphate, sulfate, or carboxylate groups, likely a biopolymer"
    if (num_amides > 4 or num_ethers > 5 or num_esters > 5) and mol_wt > 400:
        return True, "Contains multiple amides, ethers, or esters, potentially a polymer"
    
    # Check for oligomeric structures with multiple linkages
    if mol_wt > 400 and num_rotatable_bonds > 10 and (num_amides > 2 or num_ethers > 2 or num_esters > 2):
        return True, "Contains repeating or long-chain structure, potentially an oligomer"
    
    # Get molecule name by InChiKey and check if it contains 'poly'
    name = Chem.MolToInchiKey(mol)
    if name is not None:
       name = Chem.MolToInchi(mol)
       if name is not None and 'poly' in name.lower():
          return True, "Name contains poly prefix"
       name = Chem.MolToSmiles(mol)
       if name is not None and 'poly' in name.lower():
          return True, "Name contains poly prefix"
          
    return False, "Does not meet polymer criteria"