"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules arising from the oxidation of C20 EFAs.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for C20 backbone
    c20_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(c20_chain_pattern):
        return False, "No C20 backbone found."
    
    # 2. Check for Carboxylic Acid or Ester Group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxylic acid or ester group found."


    # 3. Check for relevant oxygenation patterns (hydroxyl, epoxide, hydroperoxide, ketone etc.) and/or conjugated systems

    # Check for Prostane ring + typical oxygenation: Hydroxyl or ketone groups
    prostane_pattern = Chem.MolFromSmarts("[C]1[C]([C])[C]([O])[C]([C])C1")
    prostane_ring_match = mol.GetSubstructMatches(prostane_pattern)
    if prostane_ring_match:
          # Check for prostane ring oxygenation: 
          hydroxyl_on_ring = Chem.MolFromSmarts("[C]1[C]([O])[C][C][C]1")
          ketone_on_ring = Chem.MolFromSmarts("[C]1[C](=[O])[C][C][C]1")
          if mol.HasSubstructMatch(hydroxyl_on_ring) or mol.HasSubstructMatch(ketone_on_ring):
               return True, "Matches Icosanoid criteria (C20 backbone with prostane like ring with required oxygenation)"

    # Check for Leukotriene-like features with conjugated double bonds and hydroxyl group
    leukotriene_pattern = Chem.MolFromSmarts("[C][S][C][C](O)[C]=[C][C]=[C][C]=[C][C]")  # includes a thioether and conjugated system
    if mol.HasSubstructMatch(leukotriene_pattern):
          return True, "Matches Icosanoid criteria (C20 backbone with conjugated double bonds and oxygenation)"
    
    # Check for hydroperoxide
    hydroperoxide_pattern = Chem.MolFromSmarts("[O]-[O]")
    if mol.HasSubstructMatch(hydroperoxide_pattern):
          return True, "Matches Icosanoid criteria (C20 backbone with hydroperoxide)"
    
    # Check for general conjugated systems + oxygenation (hydroxy, epoxide)
    conjugated_system_with_oxygen = Chem.MolFromSmarts("[C]=[C][C]=[C][C](O)[C]") # a conjugated system and a hydroxyl
    if mol.HasSubstructMatch(conjugated_system_with_oxygen):
          return True, "Matches Icosanoid criteria (C20 backbone with conjugated double bonds and oxygenation)"

    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Matches Icosanoid criteria (C20 backbone with epoxide)"
    
    # Check for general conjugated systems + carbonyl and hydroxyl group
    conjugated_system_with_carbonyl_and_hydroxyl = Chem.MolFromSmarts("[C]=[C][C]=[C][C](=O)[C](O)[C]") # a conjugated system and a carbonyl and hydroxyl
    if mol.HasSubstructMatch(conjugated_system_with_carbonyl_and_hydroxyl):
      return True, "Matches Icosanoid criteria (C20 backbone with conjugated double bonds and carbonyl/hydroxyl group)"
    
    
    return False, "Does not match Icosanoid criteria"