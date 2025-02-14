"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid contains the C6-C3-C6 core structure, typically as a benzopyran derivative,
    and includes various subgroups such as flavones, flavanones, isoflavones, chalcones, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for overall size, flavonoids are typically < 1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for a flavonoid"
    
    # Flexible C6-C3-C6 core: two aromatic rings linked by 3 atoms.
    # This pattern allows for heteroatoms and sp2/sp3 linkers.
    core_pattern = Chem.MolFromSmarts("([c,n,o]1[c,n,o][c,n,o][c,n,o][c,n,o][c,n,o]1)[c,C,n,o][c,C,n,o][c,C,n,o]([c,n,o]1[c,n,o][c,n,o][c,n,o][c,n,o][c,n,o]1)")
    if not mol.HasSubstructMatch(core_pattern):
      # Check specifically for chalcones (C6-C3-C6 with a carbonyl)
      chalcone_pattern = Chem.MolFromSmarts("[c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1C(=O)C[C][c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1")
      if not mol.HasSubstructMatch(chalcone_pattern):
            return False, "No C6-C3-C6 core structure found"
    
    # Check for common Flavonoid scaffolds:
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2c(cc1)occc2")  # Benzopyran
    chromene_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(=O)cc2")  # Chromene
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)cc(c2)c3ccccc3") # Isoflavone scaffold
    flavanone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)c(=O)c[cH]2c3ccccc3") #Flavanone scaffold

    if not (mol.HasSubstructMatch(benzopyran_pattern) or mol.HasSubstructMatch(chromene_pattern) or  mol.HasSubstructMatch(isoflavone_pattern) or mol.HasSubstructMatch(flavanone_pattern)):
        return False, "No common benzopyran or chromene-like ring system found"


    # Check for oxygen atoms (loosened the requirement), at least one is a must for flavonoids
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
       return False, "Not enough oxygen atoms for flavonoid."

    return True, "Likely a flavonoid based on core structure and heterocyclic ring"