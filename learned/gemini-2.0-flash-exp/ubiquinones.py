"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones have a benzoquinone core, two methoxy groups (or variations) at positions 2 and 3,
    a methyl group at position 5, and a polyprenoid side chain (or just a methyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure: benzoquinone with 2,3-dimethoxy/dihydroxy or mix, and 5-methyl
    #Accept also deprotonated oxygens. 
    core_pattern = Chem.MolFromSmarts("[c]1[c]([OX2H-])[c]([OX2H-])[c](=O)[c]([CH3])[c]1=O")
    if not mol.HasSubstructMatch(core_pattern):
         return False, "Core benzoquinone structure (2,3-dimethoxy/dihydroxy, 5-methyl) not found"

    # Check for polyprenoid side chain or methyl group at position 6. 
    #For ubiquinone-0 there is no side chain. Check for a single methyl
    side_chain_pattern = Chem.MolFromSmarts("[c]1[c]([OX2H-])[c]([OX2H-])[c](=O)[c]([CH3])[c]([CX4])1=O")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)

    no_side_chain_pattern = Chem.MolFromSmarts("[c]1[c]([OX2H-])[c]([OX2H-])[c](=O)[c]([CH3])[c]([CH3])1=O")
    no_side_chain_matches = mol.GetSubstructMatches(no_side_chain_pattern)

    if len(side_chain_matches) < 1 and len(no_side_chain_matches) < 1:
      return False, "No sidechain or methyl group in position 6 found."

    # Check if polyprenoid side chain
    branched_chain_pattern = Chem.MolFromSmarts("[CX4]([CH3])([CH3])") #Branched alkyl
    if len(side_chain_matches) > 0:
        if not mol.HasSubstructMatch(branched_chain_pattern):
             #Try to match simple methyl if no chain
           if len(no_side_chain_matches) == 0:
               return False, "Polyprenoid or methyl side chain required at pos 6"

    return True, "Matches the ubiquinone criteria: core structure with 2,3-dimethoxy/dihydroxy, 5-methyl, and polyprenoid side chain or methyl group at position 6"