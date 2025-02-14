"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    An alpha-amino acid zwitterion has a protonated amino group and a deprotonated carboxyl group on the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid zwitterion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for zwitterionic alpha-amino acid.
    # The pattern is:
    #   - [NX3;H2,H3+] (nitrogen with 2 or 3 hydrogens) - an amino group
    #   - [CX4] (a carbon with four bonds, this is the alpha carbon)
    #   - [CX3](=[OX1])[OX2-] a carboxylate group.
    # We specify direct bonds to the alpha carbon
    zwitterion_pattern = Chem.MolFromSmarts("[NX3;H2,H3+][CX4][CX3](=[OX1])[OX2-]")

    # Get substructure matches
    matches = mol.GetSubstructMatches(zwitterion_pattern)

    # Check if at least one substructure is found
    if not matches:
          return False, "Missing zwitterionic alpha-amino acid core structure"


    # Additionally check that the matching alpha carbon has ONLY an amino and carboxyl group, not some additional groups on the alpha carbon,
    # for example, a phospho group. For this, we need to get the indices of the matched atoms and check if other atoms are directly attached.
    # Loop through each match
    for match in matches:
         alpha_carbon_idx = match[1] # index of the central carbon
         alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
         
         # Get neighbors, and check that only the nitrogen and carboxyl are attached. Other neighbors must be hydrogens
         neighbors = [neighbor.GetIdx() for neighbor in alpha_carbon.GetNeighbors()]
         if len(neighbors) != 4: #if less or more than four bonds
             return False, "Too many or too few atoms bonded to the alpha carbon."

         # Check for the attached Nitrogen
         n_idx = match[0]
         if n_idx not in neighbors:
              return False, "Nitrogen not correctly attached to the alpha carbon"
              
         # Check for the attached carbon of the carboxylic group
         c_idx = match[2]
         if c_idx not in neighbors:
              return False, "Carboxylate not correctly attached to the alpha carbon"


    return True, "Molecule contains the zwitterionic alpha-amino acid structure"