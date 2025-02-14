"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.

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

    # 1. Core Phenylpropane Structure
    phenylpropane_pattern = Chem.MolFromSmarts("c1ccccc1CC[CX4]") # C-C-C, allow for any C
    if mol.HasSubstructMatch(phenylpropane_pattern):
            
        # 2. Look for typical modifications or functional groups.
        hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
        carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
        methoxy_pattern = Chem.MolFromSmarts("[CX4]-[OX2]")
        ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
        glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@H](CO)O)O[C@H]1[OX2]")
        
        # 3. Check for common heterocycles
        benzopyran_pattern = Chem.MolFromSmarts("c1ccccc1-C=2OC=CC2")
        benzodioxole_pattern = Chem.MolFromSmarts("c1ccccc1-C2-O-C-O-2")
        isobenzofuranone_pattern = Chem.MolFromSmarts("c1ccccc1-C2-C(=O)-O-2")
        coumarin_pattern = Chem.MolFromSmarts("c1ccccc1-C=2C=CC(=O)O2")
        
        #4. Check for flavone skeleton
        flavone_pattern = Chem.MolFromSmarts("c1ccccc1-c2oc(cc(=O)c2c1)-c3ccccc3")

        # Check for matches
        if (mol.HasSubstructMatch(hydroxyl_pattern) or
            mol.HasSubstructMatch(carbonyl_pattern) or
            mol.HasSubstructMatch(methoxy_pattern) or
            mol.HasSubstructMatch(ester_pattern) or
            mol.HasSubstructMatch(glycoside_pattern) or
            mol.HasSubstructMatch(benzopyran_pattern) or
            mol.HasSubstructMatch(benzodioxole_pattern) or
            mol.HasSubstructMatch(isobenzofuranone_pattern) or
            mol.HasSubstructMatch(coumarin_pattern) or
            mol.HasSubstructMatch(flavone_pattern)
        ):
              return True, "Phenylpropane core with common phenylpropanoid functional groups or heterocyclic ring system detected."
        return False, "Phenylpropane core found but no common phenylpropanoid modifications found."
    return False, "No phenylpropane core (C6-C3) structure detected."