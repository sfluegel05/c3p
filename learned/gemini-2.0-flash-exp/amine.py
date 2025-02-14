"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound derived from ammonia by replacing one, two or three
    hydrogen atoms by hydrocarbyl groups

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen atom
    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_nitrogen:
        return False, "No nitrogen atom found"

    # Check that at least one N is bonded to at least one C.
    # This handles primary, secondary, and tertiary amines, as well as quaternary
    # ammoniums.  It will also work if nitrogen has positive charge
    amine_pattern = Chem.MolFromSmarts('[NX4,NX3,NX2;!$([NX3]=[O]),!$([NX3](=[O])([!#6])!#6)]-[#6]')
    if not mol.HasSubstructMatch(amine_pattern):
         return False, "No nitrogen directly bonded to a carbon atom"

    # Exclude amides  (-C(=O)-N-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    if mol.HasSubstructMatch(amide_pattern):
      return False, "Contains amide group"

    # Exclude nitro groups (-NO2)
    nitro_pattern = Chem.MolFromSmarts("[N](=[O])=[O]")
    if mol.HasSubstructMatch(nitro_pattern):
        return False, "Contains nitro group"
        
    return True, "Contains a nitrogen atom bonded to at least one carbon atom and is not amide nor nitro."