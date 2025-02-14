"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated fatty acids that are required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for unwanted functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Contains hydroxyl group, not a typical fatty acid"
    
    peroxide_pattern = Chem.MolFromSmarts("O-O")
    if mol.HasSubstructMatch(peroxide_pattern):
        return False, "Contains peroxide group, not a typical fatty acid"

    ether_pattern = Chem.MolFromSmarts("C-O-C")
    if mol.HasSubstructMatch(ether_pattern):
         return False, "Contains ether linkage, not a typical fatty acid"

    ester_pattern = Chem.MolFromSmarts("C(=O)[OX2;!H]") #Avoid matching the carboxylic acid
    if len(mol.GetSubstructMatches(ester_pattern)) > 1: #Allow one matching carboxylic acid
      return False, "Contains ester linkage, not a typical fatty acid"
      
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide linkage, not a typical fatty acid"

    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)")
    if mol.HasSubstructMatch(phosphate_pattern):
         return False, "Contains phosphate group, not a typical fatty acid"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)")
    if mol.HasSubstructMatch(glycerol_pattern):
      return False, "Contains glycerol backbone, not a typical fatty acid"

    # Count carbons. We are more flexible now and allow for slightly longer or shorter chains.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not 14 <= c_count <= 40:
      return False, "Carbon count is not within the typical range (14-40)"
    
    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
          return False, f"Found {len(double_bond_matches)} double bonds, need at least 2"

    # Check for long aliphatic chain with at least 10 carbons. The previous version had [CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3], which was too short.
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_chain_pattern)
    if len(aliphatic_matches) < 1:
         return False, "Molecule does not contain long aliphatic chain"
    

    # Check for the presence of rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a typical fatty acid"

    # Check if the double bonds are in the chain:
    chain_atoms = []
    for match in aliphatic_matches:
      chain_atoms.extend(match)
    
    for bond_match in double_bond_matches:
      if not (mol.GetBondBetweenAtoms(bond_match[0],bond_match[1]) is not None and
          bond_match[0] in chain_atoms and bond_match[1] in chain_atoms):
              return False, "Double bond not part of main aliphatic chain"

    # Check for the presence of trans double bonds and be more strict by ensuring all double bonds are cis.
    cis_double_bond_pattern = Chem.MolFromSmarts("[CX3;!R]=[CX3;!R]")
    cis_matches = mol.GetSubstructMatches(cis_double_bond_pattern)

    if len(cis_matches) != len(double_bond_matches): #All double bonds must be cis
      return False, "All double bonds need to be cis"
    

    # Check molecular weight range (adjusting the range a bit to be more flexible)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not 200 <= mol_wt <= 700:
       return False, "Molecular weight is outside of the typical range for an essential fatty acid (200-700)"


    return True, "Contains a carboxylic acid group, a long linear carbon chain (C14-C40), multiple cis double bonds (at least 2), and a suitable molecular weight, without unwanted groups"