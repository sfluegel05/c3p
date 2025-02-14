"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    Bile acid conjugates have a steroidal core and are conjugated to a hydrophilic/charged group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general bile acid core pattern (4 fused rings)
    bile_acid_core_pattern = Chem.MolFromSmarts("C1C2C3C(C)(C1)CCC2C4(C)CCC3C4")
    if not mol.HasSubstructMatch(bile_acid_core_pattern):
         return False, "No bile acid core structure detected."
     
    # Define common conjugation groups using SMARTS patterns

    # Glycine, taurine, other aminoacids
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4]([CX4,OX2])[CX3](=[OX1])") # Very broad, needs improvement
    
    # Taurine specifically (amino sulfonic acid):
    taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)[O,H]")
    
    # Sulfates
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,H]")

    # Glucuronides
    glucuronide_pattern = Chem.MolFromSmarts("C1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)C(=O)O")
    
    # Sugars
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1")
    
    # Coenzyme A
    coenzymeA_pattern = Chem.MolFromSmarts("C(CC(=O)NCC(CC(C)O)C(O)P(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)OP(=O)(O)O)C(O)(C)C)S")
   
   # Check for the presence of at least one conjugation group
    if not (mol.HasSubstructMatch(amino_acid_pattern) or \
            mol.HasSubstructMatch(taurine_pattern) or \
            mol.HasSubstructMatch(sulfate_pattern) or \
            mol.HasSubstructMatch(glucuronide_pattern) or \
            mol.HasSubstructMatch(sugar_pattern) or \
            mol.HasSubstructMatch(coenzymeA_pattern)):
        return False, "No common conjugation group detected."
    
    # Count the number of carbons and oxygens as a sanity check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for a bile acid conjugate"
    
    if o_count < 3:
         return False, "Too few oxygens for a bile acid conjugate"


    return True, "Contains a bile acid core and at least one common conjugation group."