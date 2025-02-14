"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    Saturated fatty acids have a carboxylic acid group (-COOH) and only carbon-carbon single bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for absence of C=C or C#C bonds
    unsaturated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]") # checking for double bonds
    if mol.HasSubstructMatch(unsaturated_pattern):
          return False, "Contains carbon-carbon double bonds"
    
    triple_pattern = Chem.MolFromSmarts("[CX2]#[CX2]") # checking for triple bonds
    if mol.HasSubstructMatch(triple_pattern):
          return False, "Contains carbon-carbon triple bonds"

    # Check for at least 3 carbons in the chain using rdkit tools.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbons to be a fatty acid"

    # Check that the structure does not have more than one carboxylic acid group.
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, must be exactly 1."
    
    # Check if rings are present (and if so, only carbons)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
       for ring in ring_info.AtomRings():
            for atom_idx in ring:
                if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6:
                   return False, "Non-carbocyclic ring present"
    

    return True, "Contains a carboxylic acid group and only carbon-carbon single bonds"