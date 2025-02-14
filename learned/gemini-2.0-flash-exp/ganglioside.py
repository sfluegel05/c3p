"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid (ceramide and oligosaccharide) with one or more sialic acids linked on the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Ceramide Identification
    # Broader pattern for ceramide - captures common variations
    ceramide_pattern = Chem.MolFromSmarts("[CHX4][CHX4](O)[CHX4](-[OX2])~[NX3][CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide structure found."

    # 2. Sialic Acid (Neu5Ac and Neu5Gc) Detection
    neu5ac_pattern = Chem.MolFromSmarts("CC(=O)N[C][C](C(=O)O)([C])([C]([C]([C]([C])O)O)O)O")
    neu5gc_pattern = Chem.MolFromSmarts("C(O)N[C][C](C(=O)O)([C])([C]([C]([C]([C])O)O)O)O")

    has_sialic_acid = mol.HasSubstructMatch(neu5ac_pattern) or mol.HasSubstructMatch(neu5gc_pattern)

    if not has_sialic_acid:
        return False, "No sialic acid (Neu5Ac or Neu5Gc) found."


    # 3. Oligosaccharide Detection (looking for multiple sugar rings)
    sugar_pattern = Chem.MolFromSmarts("C1[OX2]C([CX4])C([CX4])C([CX4])C1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 2 :
        return False, "No oligosaccharide chain found."


    # Basic confirmation of link between Ceramide and Sugar Chain (check for an ether in the ceramide attached to a sugar ring)
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if ceramide_matches:
        ceramide_oxygen_idx = [match[3] for match in ceramide_matches[0]]
        if ceramide_oxygen_idx:
           ceramide_oxygen_atom = mol.GetAtomWithIdx(ceramide_oxygen_idx[0])
           linked_to_sugar = False
           for neighbor in ceramide_oxygen_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and mol.HasSubstructMatch(sugar_pattern, [neighbor.GetIdx()]):
                    linked_to_sugar = True
                    break
           if not linked_to_sugar:
             return False, "Ceramide not linked to a sugar chain"


    return True, "Contains a ceramide, an oligosaccharide, and at least one sialic acid (Neu5Ac or Neu5Gc) unit."