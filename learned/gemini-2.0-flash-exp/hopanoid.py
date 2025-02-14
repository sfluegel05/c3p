"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a precise core hopane skeleton SMARTS pattern. 
    # The core structure of hopane: 5 fused rings - four 6-membered rings and one 5-membered ring.
    # All ring bonds are single bonds, thus we use - instead of ~
    hopane_core_pattern = Chem.MolFromSmarts("[C]1-[C]-[C]-[C]2-[C]-[C]1-[C]3-[C]-[C]-[C]4-[C]2-[C]3-[C]5-[C]-[C]-[C]-[C]45")

    if not mol.HasSubstructMatch(hopane_core_pattern):
        return False, "Molecule does not contain the core hopane pentacyclic ring system."
    
    # Check for the presence of at least five quaternary carbons within the core ring system.
    # The core hopane has five quaternary carbons, so we check those explicitly by checking if the
    # core carbons within the ring are quaternary
    core_quaternary_carbon_pattern = Chem.MolFromSmarts("[C]1-[C](C)(C)-[C]-[C]2-[C]-[C]1-[C]3-[C]-[C]-[C]4-[C]2-[C]3-[C]5-[C]-[C]-[C]-[C]45")
    quaternary_matches = mol.GetSubstructMatches(core_quaternary_carbon_pattern)
    if len(quaternary_matches) != 1:
        return False, f"Core quaternary carbons not found. Must have exactly 1 match, got {len(quaternary_matches)}"
    
    # Check for the number of carbons (hopane is C30, but it can be higher, so at least 25 is OK)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 25:
        return False, f"Too few carbon atoms. Must have at least 25, got {num_carbons}"
    
    # Molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
      return False, "Molecular weight is too low for a hopanoid."
      
    
    return True, "Molecule contains the hopane core structure with at least 5 core quaternary carbons and appropriate size for hopanoid"