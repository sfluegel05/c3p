"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    #Handle mixtures and salts
    smiles_parts = smiles.split('.')
    if len(smiles_parts) > 1:
       for part in smiles_parts:
         mol = Chem.MolFromSmiles(part)
         if mol is None:
            return False, "Invalid SMILES string"
       # Mixture and salt: if one is a large molecule with a repeating pattern, it may be a polymer, otherwise return false.
       for part in smiles_parts:
         mol = Chem.MolFromSmiles(part)
         if mol is not None:
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            if mol_wt > 500 and num_rotatable_bonds > 10: #check for longer chains or cyclic structures.
                 return True, "Mixture containing large molecule with long chains, potentially a polymer" #at least one large molecule with long chain
       return False, "Mixture of small molecules, not a polymer"
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic polymer criteria
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    #Check for chains, but also allow cyclic, highly substituted molecules, and those with heteroatoms.
    if mol_wt > 500 and num_rotatable_bonds > 10 :
        return True, "High molecular weight and long chain, possibly a polymer"

    #Check for phosphate or sulfate groups. Some biopolymers contain these groups, such as nucleic acids.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)(O)")
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    
    if mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(sulfate_pattern) or mol.HasSubstructMatch(carboxylate_pattern):
        if mol_wt > 300: #some biopolymers can be smaller, but many are larger
            return True, "Contains phosphate, sulfate or carboxylate groups, and is a larger molecule, potentially a biopolymer"

    
    # Check for repeating units or long chains:
    chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]") # 10 carbon chain
    if mol.HasSubstructMatch(chain_pattern):
        return True, "Contains long carbon chain, possibly a polymer"

    #Check for explicit "poly" in the name.
    
    name = Chem.MolToInchiKey(mol)
    if name is not None and 'poly' in name.lower():
       return True, "Name contains poly prefix"

    return False, "Does not meet polymer criteria"