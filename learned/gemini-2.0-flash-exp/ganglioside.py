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
    ceramide_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CHX4](-[OX2])~[CHX3]~[NX3]~[CX3](=[OX1])-[CX4]")
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide structure found."

    # 2. Sialic Acid (Neu5Ac and Neu5Gc) Detection
    neu5ac_pattern = Chem.MolFromSmarts("[C]([C](=[O])N[C]1[C]([CH]([CH]([CH]([C]([CH]([CH]([CH]1O)O)O)O)C(=O)O)O)O)O")
    neu5gc_pattern = Chem.MolFromSmarts("[C]([C](=[O])N[C]1[C]([CH]([CH]([CH]([C]([CH]([CH]([CH]1O)O)O)O)C(=O)O)O)O)CO")
    neu5ac_matches = mol.GetSubstructMatches(neu5ac_pattern)
    neu5gc_matches = mol.GetSubstructMatches(neu5gc_pattern)

    if not neu5ac_matches and not neu5gc_matches:
        return False, "No sialic acid (Neu5Ac or Neu5Gc) found."

    # 3. Check for Glycosidic bonds (ether oxygens) - ensure sugar is present
    ether_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and len(atom.GetBonds()) > 1 )
    if ether_count < 3:
       return False, "No sugar chain found."

    # 4. Verify ceramide link to sugar and sialic acid link to sugar
    # Check if the -OH in the ceramide is linked to a sugar, which has an ether oxygen.
    ceramide_atom_idx = [match[2] for match in ceramide_matches[0]] # index of the oxygen atom in the ceramide
    if ceramide_atom_idx:
        ceramide_atom = mol.GetAtomWithIdx(ceramide_atom_idx[0])
        has_sugar_link = False
        for neighbor in ceramide_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and len(neighbor.GetBonds()) > 1:
                has_sugar_link = True
                break
        if not has_sugar_link:
            return False, "Ceramide not linked to a sugar."

    # 5. Check for multiple pyranose rings as a sign of oligosaccharide
    pyranose_pattern = Chem.MolFromSmarts("C1OC(C)C(C)C(C)C1")
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    if len(pyranose_matches) < 2:
        return False, "Not an oligosaccharide chain."


    return True, "Contains a ceramide, an oligosaccharide, and at least one sialic acid (Neu5Ac or Neu5Gc) unit."