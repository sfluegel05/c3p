"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids have a phenylpropane skeleton (C6-C3 or C6-C2) and various modifications

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core Phenylpropane Structure (C6-C3 or C6-C2)
    # More specific patterns to prevent matching other random carbon chains on phenyl
    phenylpropane_C3_pattern = Chem.MolFromSmarts("c1ccccc1[CH2][CH2][CH2]")
    phenylpropane_C2_pattern = Chem.MolFromSmarts("c1ccccc1[CH2][CH2]")
    phenylpropane_C2_O_pattern = Chem.MolFromSmarts("c1ccccc1[CH2][O][CH2]") # more specific pattern to prevent false matches


    # 2. Common Modifications (hydroxyl, carbonyl, methoxy, ester etc)
    modification_pattern1 = Chem.MolFromSmarts("[OH]")
    modification_pattern2 = Chem.MolFromSmarts("[CX3]=[OX1]")
    modification_pattern3 = Chem.MolFromSmarts("[CX4]-[OX2]")
    modification_pattern4 = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    modification_pattern5 = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@H](CO)O)O[C@H]1[OX2]") # glycoside

    # 3. Common phenylpropanoid ring systems:
    # More generic SMARTS patterns
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc1-C=2OC=CC2") # more generic
    benzodioxole_pattern = Chem.MolFromSmarts("c1ccccc1-C2-O-C-O-2")
    isobenzofuranone_pattern = Chem.MolFromSmarts("c1ccccc1-C2-C(=O)-O-2")
    coumarin_pattern = Chem.MolFromSmarts("c1ccccc1-C=2C=CC(=O)O2") # more general pattern to match different coumarins
    flavone_pattern = Chem.MolFromSmarts("c1ccccc1-c2oc(cc(=O)c2)-c3ccccc3") # more general pattern
    flavonoid_core_pattern = Chem.MolFromSmarts("c1ccccc1c2c(=O)cc(O)c(-c3ccccc3)o2")
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)c2ccccc2")
    dihydrobenzofuran_pattern = Chem.MolFromSmarts("c1ccccc1[CH2]2OC[CH2]2") # pattern for isograndidentatin type molecules
    isoflavone_pattern = Chem.MolFromSmarts("c1ccccc1c2c(=O)cc(O)c(c3ccccc3)o2") #isoflavone

    # Check for core and modification (AND condition).
    # Prioritize checking ring systems.
    if (benzopyran_pattern and mol.HasSubstructMatch(benzopyran_pattern)) or \
            (benzodioxole_pattern and mol.HasSubstructMatch(benzodioxole_pattern)) or \
            (isobenzofuranone_pattern and mol.HasSubstructMatch(isobenzofuranone_pattern)) or \
            (coumarin_pattern and mol.HasSubstructMatch(coumarin_pattern)) or \
             (flavone_pattern and mol.HasSubstructMatch(flavone_pattern)) or \
            (flavonoid_core_pattern and mol.HasSubstructMatch(flavonoid_core_pattern)) or \
            (chalcone_pattern and mol.HasSubstructMatch(chalcone_pattern)) or \
            (dihydrobenzofuran_pattern and mol.HasSubstructMatch(dihydrobenzofuran_pattern)) or \
            (isoflavone_pattern and mol.HasSubstructMatch(isoflavone_pattern)):
         return True, "Phenylpropanoid ring system detected."

    has_core = False
    if mol.HasSubstructMatch(phenylpropane_C3_pattern) or \
       mol.HasSubstructMatch(phenylpropane_C2_pattern) or \
       mol.HasSubstructMatch(phenylpropane_C2_O_pattern):
        has_core = True

    if has_core:
         has_modification = False
         if (modification_pattern1 and mol.HasSubstructMatch(modification_pattern1)) or \
           (modification_pattern2 and mol.HasSubstructMatch(modification_pattern2)) or \
            (modification_pattern3 and mol.HasSubstructMatch(modification_pattern3)) or \
            (modification_pattern4 and mol.HasSubstructMatch(modification_pattern4)) or \
            (modification_pattern5 and mol.HasSubstructMatch(modification_pattern5)):
              has_modification = True
         if has_modification:
            return True, "Phenylpropane core with typical modification detected."
         else:
            return True, "Phenylpropane core detected, no specific modification found but core detected."
    
    return False, "No phenylpropane core (C6-C2 or C6-C3) structure detected."